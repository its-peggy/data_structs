#include <vector>
#include <iostream>

using namespace std;

template<class T>
struct Link {
    explicit Link(const T& info, Link *next = 0) : info(info), next(next) { }
    // This avoids stack overflow
    ~Link() {
        Link *p = next;
        while (p) {
            Link *q = p->next;
            p->next = 0;
            delete p;
            p = q;
        }
    }
    T info;
    Link *next;
};

vector<int> josephus(int n, int k) {
    vector<int> v;
    int move;
    // create circular LL
    Link<int>* tail = new Link<int>(n);
    Link<int>* head = tail;
    for (int i = n-1; i >= 1; --i) {
        Link<int>* temp = new Link<int>(i, head);
        head = temp;
    }
    tail->next = head;
    // traverse
    Link<int>* cur = tail;
    int left = n;
    while (cur->next != cur) {
        move = (k-1) % left;
        for (int i = 1; i <= move; ++i) {
            cur = cur->next;
        }
        Link<int>* kill = cur->next;
        cur->next = kill->next;
        v.push_back(kill->info);
        kill->next = nullptr;
        delete kill;
        left--;
    }
    v.push_back(cur->info);
    return v;
}

template<class T>
vector<int> loopTail(Link<T>* head) {
    vector<int> v;
    // check if linked list is empty
    if (!head) {
        v.push_back(0);
        v.push_back(0);
        return v;
    }
    // STAGE 1
    bool cycle = true;
    Link<T>* slow = head;
    Link<T>* fast = head;
    do {
        if (slow->next != nullptr) { slow = slow->next; }
        else { cycle = false; break; }
        if (fast->next != nullptr) { fast = fast->next; }
        else { cycle = false; break; }
        if (fast->next != nullptr) { fast = fast->next; }
        else { cycle = false; break; }
    }
    while (slow != fast);

    // if no cycle
    if (!cycle) {
        int tail_len = 0;
        Link<T>* start = head;
        while (start->next != nullptr) {
            start = start->next;
            tail_len++;
        }
        v.push_back(0);
        v.push_back(tail_len);
        return v;
    }
    // STAGE 2 -- exist cycle, slow=fast at this point
    int loop_len = 0;
    do {
        fast = fast->next;
        loop_len++;
    }
    while (slow != fast);

    v.push_back(loop_len);

    // STAGE 3
    slow = head;
    fast = head;
    for (int i = 0; i < loop_len; ++i) {
        fast = fast->next;
    }

    int tail_len = 0;
    while (slow != fast) {
        slow = slow->next;
        fast = fast->next;
        tail_len++;
    }
    v.push_back(tail_len);
    return v;
}
