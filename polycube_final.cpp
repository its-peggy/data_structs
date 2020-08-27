/*
Feb 2020
Generate solutions to 3D polycube puzzles
Ref: Knuth's "Dancing Links"
https://arxiv.org/pdf/cs/0011047.pdf
*/

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <utility>
#include <stack>
#include <string>
#include <cctype>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <initializer_list>

using namespace std;
using namespace std::chrono;

/////////////////////////////////////// structs + typedefs /////////////////////////////////////////

struct Point {
    int x;
    int y;
    int z;
    Point(int x1, int y1, int z1) {
        x = x1;
        y = y1;
        z = z1;
    }
    Point(const Point &p2) {
        x = p2.x;
        y = p2.y;
        z = p2.z;
    }
};

struct Column {
    Column* up;
    Column* down;
    Column* left;
    Column* right;
    Column* column;
    int size;
    string name;
    Column() {
        up = nullptr;
        down = nullptr;
        left = nullptr;
        right = nullptr;
        column = nullptr;
        size = 0;
        name = "";
    }
    Column(string s) {
        up = nullptr;
        down = nullptr;
        left = nullptr;
        right = nullptr;
        column = nullptr;
        size = 0;
        name = s;
    }
};

typedef struct Column Column;

typedef vector<Point> Piece;
typedef vector<vector<Point> > VPiece;
typedef vector<vector<vector<Point> > > VVPiece;
typedef vector<Column*> VC;
typedef vector<vector<Column*> > VVC;

typedef vector<vector<int> > VVI;
typedef vector<vector<vector<int> > > VVVI;

//////////////////////////////////////// global variables //////////////////////////////////////////

int N, D, W, H;
vector<vector<int> > preM;
VVC M;
Column* root;
VC solutions;
vector<string> final_solutions;
unordered_set<string> known, distinct, real_ans;

string rotations_arr[] = { "X", "Y", "Z", "XX", "XY", "XZ", "YX", "YY", "ZY", "ZZ", "XXX", "XXY", "XXZ", "XYX", "XYY", "XZZ", "YXX", "YYY", "ZZZ", "XXXY", "XXYX", "XYXX", "XYYY" };
vector<string> rotations(begin(rotations_arr), end(rotations_arr));

string perms_arr[] = { "012", "021", "102", "120", "201", "210" };
vector<string> perms(begin(perms_arr), end(perms_arr));

string negs_arr[] = { "000", "001", "010", "100", "011", "101", "110", "111" };
vector<string> negs(begin(negs_arr), end(negs_arr));

///////////////////////////////////////// << overloading ///////////////////////////////////////////

ostream& operator<<(ostream& os, const Point& p) {
    os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

/////////////////////////////////////// rotation functions /////////////////////////////////////////

// (2, 1, 0) --> (-1, 2, 0)
void rot90cc_z(Piece &piece) {
    for (int i = 0; i < piece.size(); ++i) {
        int temp_x = piece[i].x;
        int temp_y = piece[i].y;
        piece[i].x = temp_y;
        piece[i].y = (-1) * temp_x;
    }
}

// (-1, 0, 2) --> (-2, 0, -1)
void rot90cc_y(Piece &piece) {
    for (int i = 0; i < piece.size(); ++i) {
        int temp_x = piece[i].x;
        int temp_z = piece[i].z;
        piece[i].x = (-1) * temp_z;
        piece[i].z = temp_x;
    }
}

// (0, 2, 1) --> (0, -1, 2)
void rot90cc_x(Piece &piece) {
    for (int i = 0; i < piece.size(); ++i) {
        int temp_y = piece[i].y;
        int temp_z = piece[i].z;
        piece[i].y = (-1) * temp_z;
        piece[i].z = temp_y;
    }
}

//////////////////////////// setup functions to generate pre-matrix ////////////////////////////////

// input: a piece
// result: the piece is moved to the 1st octant, aligned with origin
void translate_to_origin(Piece &piece) {
    int min_x = piece[0].x;
    int min_y = piece[0].y;
    int min_z = piece[0].z;
    for (int i = 1; i < piece.size(); ++i) {
        if (piece[i].x < min_x) min_x = piece[i].x;
        if (piece[i].y < min_y) min_y = piece[i].y;
        if (piece[i].z < min_z) min_z = piece[i].z;
    }
    for (int i = 0; i < piece.size(); ++i) {
        piece[i].x += (-1) * min_x;
        piece[i].y += (-1) * min_y;
        piece[i].z += (-1) * min_z;
    }
}

// input: ONE piece that has been already translated to origin
// output: Point representing max coords of any point in piece
bool in_bounds(Piece &piece) {
    int max_x = piece[0].x;
    int max_y = piece[0].y;
    int max_z = piece[0].z;
    for (int i = 0; i < piece.size(); ++i) {
        if (piece[i].x > max_x) max_x = piece[i].x;
        if (piece[i].y > max_y) max_y = piece[i].y;
        if (piece[i].z > max_z) max_z = piece[i].z;
    }
    return (max_x < D && max_y < W && max_z < H);
}

// input: vector of all 24 rotations of ONE piece, already at origin (vector of vectors of Points)
// result: duplicates are removed from configs
void remove_duplicates(VPiece &configs) {
    // one by one
    // convert to encoded coords
    // sort
    // put in unordered_set<vector<int> >
    // clear configs
    // decode each element in set
    // put decoded element into configs
    set<vector<int> > s; // each vector<int> is ONE PIECE
    for (int i = 0; i < configs.size(); ++i) { // for each rotated piece
        vector<int> encoded_coords;
        for (int j = 0; j < configs[i].size(); ++j) { // for each point in piece
            // ENCODE
            int encoding = (configs[i][j].x)*(W*H) + (configs[i][j].y)*(H) + (configs[i][j].z);
            encoded_coords.push_back(encoding);
        }
        sort(encoded_coords.begin(), encoded_coords.end());
        s.insert(encoded_coords);
    }
    configs.clear();
    // DECODE
    for (const vector<int> &vi : s) {
        Piece piece;
        for (int i = 0; i < vi.size(); ++i) {
            // D, W, H
            int x = vi[i]/(int)(W*H);
            int y = (vi[i]%(W*H))/(int)H;
            int z = vi[i]%H;
            Point p(x, y, z);
            piece.push_back(p);
        }
        configs.push_back(piece);
    }
}

// input: ONE origin-aligned piece
// result: populates all_rotations (currently empty) w/ all DISTINCT rotations.
void generate_rotations(VPiece &all_rotations, Piece &piece) {
    translate_to_origin(piece);
    if (in_bounds(piece)) {
        all_rotations.push_back(piece);
    }
    for (int i = 0; i < rotations.size(); ++i) { // generate all 24 rotations
        string s = rotations[i];
        Piece temp = piece; // start with original
        for (int j = s.length()-1; j >=0 ; --j) {
            if (s[j] == 'X') rot90cc_x(temp);
            if (s[j] == 'Y') rot90cc_y(temp);
            if (s[j] == 'Z') rot90cc_z(temp);
        }
        translate_to_origin(temp);
        if (in_bounds(temp)) {
            all_rotations.push_back(temp);
        }
    }
    remove_duplicates(all_rotations);
}

// W, D, H
// when this is called, output contains only the DISTINCT rotations of ONE piece
void generate_translations(VPiece &output) {
    int num_of_rots = output.size();
    for (int i = 0; i < num_of_rots; ++i) { // generate translations for EACH ROTATED PIECE
        VPiece all_translations; // of this piece
        Piece piece = output[i];
        int max_x = 0;
        int max_y = 0;
        int max_z = 0;
        for (int j = 0; j < piece.size() ; ++j) {
            if (piece[j].x > max_x) max_x = piece[j].x;
            if (piece[j].y > max_y) max_y = piece[j].y;
            if (piece[j].z > max_z) max_z = piece[j].z;
        }
        for (int dx = 0; dx < D-max_x; ++dx) {
            for (int dy = 0; dy < W-max_y; ++dy) {
                for (int dz = 0; dz < H-max_z; ++dz) {
                    if (dx!=0 || dy!=0 || dz!=0) {
                        Piece translated;
                        for (int j = 0; j < piece.size(); ++j) {
                            Point new_pt(piece[j].x+dx, piece[j].y+dy, piece[j].z+dz);
                            translated.push_back(new_pt);
                        }
                        all_translations.push_back(translated);
                    }
                }
            }
        }
        output.insert(output.end(), all_translations.begin(), all_translations.end());
    }
}

// input: vector of Pieces
// output: vector of (vectors of ints). each (vector of ints) is one piece.
VVI encode_VPiece(VPiece &vp) {
    VVI res;
    for (int i = 0; i < vp.size(); ++i) { // for each piece
        vector<int> cur_piece;
        for (int j = 0; j < vp[i].size(); ++j) { // for each point in current piece
            int encoding = (vp[i][j].x)*(W*H) + (vp[i][j].y)*(H) + (vp[i][j].z);
            cur_piece.push_back(encoding);
            // do not have to sort
        }
        res.push_back(cur_piece);
    }
    return res;
}

// input: ONE input string representing ONE piece
// result: adds to a vector containing all DISTINCT rotations/translations of piece
void parse(VPiece &output, string s) {
    int min_x, min_y, min_z;
    // construct piece
    Point cur(0, 0, 0);
    vector<Point> path;
    Piece piece;
    path.push_back(cur);
    piece.push_back(cur);
    for (int i = 0; i < s.length(); ++i) {
        char c = s[i];
        if (!isdigit(c)) {
            if (c == 'E') cur.x++;
            else if (c == 'W') cur.x--;
            else if (c == 'N') cur.y++;
            else if (c == 'S') cur.y--;
            else if (c == 'F') cur.z++;
            else if (c == 'B') cur.z--;
            if (cur.x < min_x) min_x = cur.x;
            if (cur.y < min_y) min_y = cur.y;
            if (cur.z < min_z) min_z = cur.z;
            path.push_back(cur);
            piece.push_back(cur);
        }
        else {
            int og_size = path.size();
            int back = c - '0';
            for (int i = og_size-1; i >= og_size-back; --i) {
                Point temp = path[i-1];
                path.push_back(temp);
                cur = temp;
            }
        }
    }
    // this function populates output (currently empty)
    generate_rotations(output, piece);
}

// input: all configs (rotations AND translations) of ONE piece, and index of cur piece
// result: matrix preM is appended with rows that represent these configs, encoded
void build_preM(VVI &configs, int index) {
    // for every piece in this group, the prefix is the same
    vector<int> prefix(N, 0);
    prefix[index] = 1;
    // the postfix depends on the actual coords
    for (int i = 0; i < configs.size(); ++i) { // loop over all rotations/trans of this ONE piece
        vector<int> pre = prefix;
        vector<int> post(D*W*H, 0);
        // populate post with 1's where they should be
        for (int j = 0; j < configs[i].size(); ++j) {
            int temp = configs[i][j];
            post[temp]++;
        }
        // puts encoded part onto back of row
        pre.insert(pre.end(), post.begin(), post.end());
        // add row to preM
        preM.push_back(pre);
    }
}

/////////////////////////////////////// building matrix M //////////////////////////////////////////

void buildM() {
    VC header_row;
    // header row = vector of column objects
    // making header row -- pt. 1
    char ch = 'A';
    string s(1, ch);
    Column* first = new Column(s);
    Column* cur = first;
    Column* prev = first;
    first->up = first;
    first->down = first;
    first->column = first;
    header_row.push_back(first);
    for (int i = 1; i < N; ++i) {
        ch++;
        string s(1, ch); // cast char ch to string
        prev = cur;
        cur = new Column(s);
        header_row.push_back(cur);
        cur->left = prev;
        prev->right = cur;
        cur->up = cur;
        cur->down = cur;
        cur->column = cur;
    }
    // making header row -- pt. 2
    for (int i = 0; i < D*W*H ; ++i) {
        string s = to_string(i);
        prev = cur;
        cur = new Column(s);
        header_row.push_back(cur);
        cur->left = prev;
        prev->right = cur;
        cur->up = cur;
        cur->down = cur;
        cur->column = cur;
    }
    // making "root" column object
    string r = "root";
    root = new Column(r);
    cur->right = root;
    root->left = cur;
    root->right = first;
    first->left = root;
    header_row.push_back(root);
    // add header row to matrix
    M.push_back(header_row);
    // making other rows
    for (int i = 0; i < preM.size(); ++i) { // iterate over each row
        string cur_row_name;
        VC new_row(N+D*W*H, nullptr);
        Column* prev = nullptr;
        for (int j = 0; j < preM[i].size(); ++j) { // iterate over this row in preM
            if (preM[i][j] == 0) { // nothing here, nullptr
                continue;
            }
            else {
                Column* thing = new Column();
                if (!prev) { // first one in row
                    cur_row_name = M[0][j]->name;
                    thing->name = cur_row_name;
                    thing->left = thing;
                    thing->right = thing;
                    thing->column = M[0][j];
                    thing->down = M[0][j];
                    thing->up = thing->down->up;
                    thing->down->up = thing;
                    thing->up->down = thing;
                }
                else { // not first in row
                    thing->name = cur_row_name;
                    thing->left = prev;
                    thing->right = prev->right;
                    prev->right = thing;
                    thing->column = M[0][j];
                    thing->down = M[0][j];
                    thing->up = thing->down->up;
                    thing->down->up = thing;
                    thing->up->down = thing;
                    thing->right->left = thing;
                    thing->left->right = thing;
                }
                thing->column->size++;
                prev = thing;
                new_row[j] = thing;
            }
            M.push_back(new_row);
        }
    }
}

////////////////////////////////////////////// DLX /////////////////////////////////////////////////

string print_solution() {
    vector<string> sol(D*W*H, "");
    Column* ok = solutions[1];
    Column* start = ok;
    while (ok->right != start) {
        ok = ok->right;
    }
    ok = ok->right;
    for (int i = 0; i < solutions.size(); ++i) { // collectively, this populates sol
        Column* cur = solutions[i];
        Column* start_right = cur;
        do {
            string s = cur->column->name;
            char n = s[0];
            if (isalpha(n)) { // is a "prefix" column -- not needed for solution string
                cur = cur->right;
            }
            else {
                string coords = cur->column->name; // 0 to D*W*H-1
                int coords_int = atoi(coords.c_str()); // cast coords to int
                sol[coords_int] = cur->name; // "A", "B", etc.
                cur = cur->right;
            }
        }
        while (cur != start_right);
    }
    string final = "";
    for (int i = 0; i < sol.size(); ++i) {
        final += sol[i];
    }
    return final;
}

Column* choose_column() {
    Column* cur = root->right; // first
    Column* res = cur;
    int cur_min = cur->size;
    while (cur->right != root) {
        if (cur->right->size < cur_min) {
            cur_min = cur->right->size;
            res = cur->right;
        }
        cur = cur->right;
    }
    return res;
}

void cover(Column*& c) {
    c->right->left = c->left;
    c->left->right = c->right;
    Column* cur1;
    Column* cur2;
    Column* start_down = c;
    while (c->down != start_down) {
        cur1 = c->down;
        Column* start_right = cur1;
        while (cur1->right != start_right) {
            cur2 = cur1->right;
            cur2->down->up = cur2->up;
            cur2->up->down = cur2->down;
            cur2->column->size -= 1;
            cur1 = cur1->right;
        }
        cur1 = cur1->right; // end at start_right
        c = c->down;
    }
    c = c->down; // end at start_down
}

void uncover(Column*& c) {
    Column* cur1;
    Column* cur2;
    Column* start_up = c;
    while (c->up != start_up) {
        cur1 = c->up;
        Column* start_left = cur1;
        while (cur1->left != start_left) {
            cur2 = cur1->left;
            cur2->column->size += 1;
            cur2->down->up = cur2;
            cur2->up->down = cur2;
            cur1 = cur1->left;
        }
        cur1 = cur1->left; // end at start_left
        c = c->up;
    }
    c = c->up;  // end at start_up
    c->right->left = c;
    c->left->right = c;
}

void search(int k) {
    if (root->right == root) {
        string s = print_solution();
        final_solutions.push_back(s);
        return;
    }
    Column* c = choose_column();
    if(c->size == 0) return;
    cover(c);
    Column* cur;
    Column* start_down = c;
    while (c->down != start_down) {
        cur = c->down;
        solutions[k] = cur;
        Column* start_right = cur;
        while (cur->right != start_right) {
            cover(cur->right->column);
            cur = cur->right;
        }
        search(k+1);
        cur = solutions[k];
        Column* start_left = cur;
        while (cur->left != start_left) {
            uncover(cur->left->column);
            cur = cur->left;
        }
        c = c->down;
    }
    c = start_down;
    uncover(c);
    return;
}

///////////////////////////// remove solution rotations/reflections ////////////////////////////////

string hash_by_order(string s) {
    unordered_map<char, char> umap;
    vector<char> res(s.length(), '\0');
    char cur = 'A';
    for (int i = 0; i < s.length(); ++i) {
        char c = s[i];
        if (umap.count(c) == 0) {
            res[i] = cur;
            umap[c] = cur;
            cur++;
        }
        else {
            res[i] = umap.at(c);
        }
    }
    string res_str(res.begin(), res.end());
    return res_str;
}

// add all 48 orientations to KNOWN
void add_rots_refs(string s) {
    // start: D x W x H
    vector<int> start(3, 0);
    start[0] = D;
    start[1] = W;
    start[2] = H;
    for (int i = 0; i < perms.size(); ++i) {
        string cur_p = perms[i];
        vector<int> cur_perm(3, 0);
        cur_perm[0] = cur_p[0] - '0';
        cur_perm[1] = cur_p[1] - '0';
        cur_perm[2] = cur_p[2] - '0';
        // PERMUTE (D, W, H) -> (new_D, new_W, new_H);
        int new_D = start[cur_perm[0]];
        int new_W = start[cur_perm[1]];
        int new_H = start[cur_perm[2]];
        for (int j = 0; j < negs.size(); ++j) {
            vector<char> res(s.length(), 'q');
            string cur_n = negs[j];
            vector<int> cur_neg(3, 0);
            cur_neg[0] = cur_n[0] - '0';
            cur_neg[1] = cur_n[1] - '0';
            cur_neg[2] = cur_n[2] - '0';
            // map every point in s
            for (int d = 0; d < D; ++d) {
                for (int w = 0; w < W; ++w) {
                    for (int h = 0; h < H; ++h) {
                        vector<int> start_pt(3, 0);
                        start_pt[0] = d;
                        start_pt[1] = w;
                        start_pt[2] = h;
                        // permute
                        int new_d = start_pt[cur_perm[0]];
                        int new_w = start_pt[cur_perm[1]];
                        int new_h = start_pt[cur_perm[2]];
                        // negate + translate
                        if (cur_neg[0] == 1) { // negative
                            new_d = (-1) * new_d;
                            new_d += new_D - 1;
                        }
                        if (cur_neg[1] == 1) {
                            new_w = (-1) * new_w;
                            new_w += new_W - 1;
                        }
                        if (cur_neg[2] == 1) {
                            new_h = (-1) * new_h;
                            new_h += new_H - 1;
                        }
                        int old_encoding = d*(W*H) + w*(H) + h;
                        int new_encoding = new_d*(new_W*new_H) + new_w*(new_H) + new_h;
                        res[new_encoding] = s[old_encoding];
                    }
                }
            }
            string res_str = "";
            for (int asdf = 0; asdf < res.size(); ++asdf) {
                res_str += res[asdf];
            }
            string res_str_hashed = hash_by_order(res_str);
            known.insert(res_str_hashed);
        }
    }
}

void process(string s) {
    string hashed = hash_by_order(s);
    if (known.count(hashed) == 0) {
        known.insert(hashed);
        distinct.insert(hashed);
        real_ans.insert(s);
    }
    add_rots_refs(hashed);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    VVPiece all_configs;
    vector<string> inputs;
    cin >> N;
    for (int i = 0; i < N; ++i) {
        solutions.push_back(nullptr);
    }
    for (int i = 0; i < N; ++i) {
        string s;
        cin >> s;
        inputs.push_back(s);
    }
    cin >> D >> W >> H;
    for (int i = 0; i < N; ++i) {
        VPiece configs_of_cur_piece;
        parse(configs_of_cur_piece, inputs[i]);
        all_configs.push_back(configs_of_cur_piece);
    }
    // put translations of each piece into its own vector in all_configs
    // appends all translations at the end of configs_of_cur_piece
    for (int i = 0; i < all_configs.size(); ++i) {
        generate_translations(all_configs.at(i)); // all_configs[i] = rotations of piece i
        VVI final = encode_VPiece(all_configs.at(i));
        build_preM(final, i);
    }
    buildM();
    search(0);
    for (int i = 0; i < final_solutions.size(); ++i) {
        process(final_solutions[i]);
    }
    cout << real_ans.size() << "\n";
    for (string s : real_ans) {
        cout << s << "\n";
    }
}
