/*
Nov 2019
Implement an efficient comparison and sort of encoded five-card poker hands
*/

#include <vector>
#include <iostream>
#include <cmath>
#include <utility>
#include <chrono>

using namespace std;
using namespace std::chrono;

typedef vector<pair<int,int> > VPII;

void insertion_sort_hand(VPII& v) {
    for (int i = 1; i < v.size(); ++i) {
        int j = i;
        while (j > 0 && ((v[j-1].first > v[j].first) ||
            ((v[j-1].first == v[j].first) && (v[j-1].second > v[j].second)))) {
            swap(v[j-1], v[j]);
            j--;
        }
    }
}

int compute_type(VPII& v, bool b) {
    int sz = v.size();
    if (b && sz == 5 &&
        (v[4].second-v[0].second == 4 ||
            (v[4].second == 12 && v[3].second == 3))) {
                if (v[4].second == 12 && v[3].second == 3) {
                    swap(v[4], v[3]);
                }
                return 9;                                           // straight flush
    }
    else if (sz == 2 && v[1].first == 4) {
        return 8;                                                   // four of a kind
    }
    else if (sz == 2 && v[1].first == 3) {
        return 7;                                                   // full house
    }
    else if (b) {
        return 6;                                                   // flush
    }
    else if (sz == 5 &&
        (v[4].second-v[0].second == 4 ||
            (v[4].second == 12 && v[3].second == 3))) {
                if (v[4].second == 12 && v[3].second == 3) {
                    swap(v[4], v[3]);
                }
                return 5;                                           // straight
    }
    else if (sz == 3 && v[2].first == 3) {
        return 4;                                                   // three of a kind
    }
    else if (sz == 3 && v[2].first == 2) {
        return 3;                                                   // two pairs
    }
    else if (sz == 4) {
        return 2;                                                   // pair
    }
    else if (sz == 5) {
        return 1;                                                   // high card
    }
    else return -1;                                                 // ERROR
}

int pow_13(int exp) {
    int res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= 13;
    }
    return res;
}

int decode_then_encode(int hand) {
    VPII v, new_v;
    bool all_same_suit = true;
    int prev_suit;
    for (int i = 0; i < 13; ++i) {
        pair<int, int> p (0, i);
        v.push_back(p);
    }
    for (int i = 0; i < 5; ++i) {
        if (all_same_suit && (i > 0) && (prev_suit != hand%4)) {
            all_same_suit = false;
        }
        prev_suit = hand%4;
        v[(hand%52)/4].first += 1;
        hand /= 52;
    }
    for (pair<int,int> p : v) {
        if (p.first != 0) new_v.push_back(p);
    }
    insertion_sort_hand(new_v);
    int type = compute_type(new_v, all_same_suit);
    int res = type * pow_13(5);;
    for (int i = 0; i < new_v.size(); ++i) {
        res += new_v[i].first * (new_v[i].second * pow_13(i));
    }
    return res;
}

//////////////////////////////////////////// SORTING /////////////////////////////////////////////

int partition(vector<int>& v, int piv, int low, int high, vector<int>& aux) {
    int i = low-1;
    int j = high+1;
    while (low <= high) {
        do {
            i++;
        } while(v[i] < piv);
        do {
            j--;
        } while(v[j] > piv);
        if (i >= j) {
            return j;
        }
        swap(v[i], v[j]);
        swap(aux[i], aux[j]);
    }
}

int median_of_three(vector<int>& v, int l, int h) {
    int a = v[l];
    int b = v[(l+h)/2];
    int c = v[h];
    if ((a <= b && b <= c) || (a >= b && b >= c)) {
        return b;
    }
    else if ((a <= c && c <= b) || (a >= c && c >= b)) {
        return c;
    }
    else return a;
}

void insertionsort(vector<int>& v, vector<int>& aux) {
    for (int i = 1; i < v.size(); ++i) {
        int j = i;
        while (j > 0 && v[j-1] > v[j]) {
            swap(v[j-1], v[j]);
            swap(aux[j-1], aux[j]);
            j--;
        }
    }
}

void quicksort(vector<int>& v, int low, int high, vector<int>& aux) {
    if (low < high) {
        if (v.size() < 30) {
            insertionsort(v, aux);
            return;
        }
        int piv = median_of_three(v, low, high);
        int p = partition(v, piv, low, high, aux);
        quicksort(v, low, p, aux);
        quicksort(v, p+1, high, aux);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void poker_sort(vector<int>& a) {
    vector<int> final_encoded;
    for (int i : a) {
        final_encoded.push_back(decode_then_encode(i));
    }
    quicksort(final_encoded, 0, final_encoded.size()-1, a);
}
