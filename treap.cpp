#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <climits>
#include <random>
#include <memory>
#include <functional>

using namespace std;

template<class T>
struct Node {
    T key;
    unsigned long int priority;
    Node* left;
    Node* right;

    Node(T key) {
        this->key = key;
        mt19937 mt_rand(time(0));
        this->priority = mt_rand()*2;
        this->left = nullptr;
        this->right = nullptr;
    }
    ~Node();
};

template<class T>
void right_rotation(Node<T>*& pivot) {
    Node<T>* l = pivot->left;
    Node<T>* r = pivot->right;
    Node<T>* ll = l->left;
    Node<T>* lr = l->right;
    pivot->left = lr;
    l->right = pivot;
    pivot = l;
}

template<class T>
void left_rotation(Node<T>*& pivot) {
    Node<T>* l = pivot->left;
    Node<T>* r = pivot->right;
    Node<T>* rl = r->left;
    Node<T>* rr = r->right;
    pivot->right = rl;
    r->left = pivot;
    pivot = r;
}

template<class T>
Node<T>* find(Node<T>*& root, T searchkey) {
    if (!root || root->key == searchkey) {
        return root;
    }
    if (root->key < searchkey) {
        return find(root->right, searchkey);
    }
    if (root->key > searchkey) {
        return find(root->left, searchkey);
    }
}

template<class T>
void insert(Node<T>*& root, T key) {
    // if empty tree
    if (!root) {
        root = new Node<T>(key);
        return;
    }
    // regular BST insertion
    if (key <= root->key) {
        insert(root->left, key);
    }
    else {
        insert(root->right, key);
    }
    // rotate to maintain heap order
    if (root->left && (root->priority < root->left->priority)) {
        right_rotation(root);
    }
    if (root->right && (root->priority < root->right->priority)) {
        left_rotation(root);
    }
}

template<class T>
Node<T>* remove(Node<T>*& root, T deletekey) { // returns root of tree after deletion
    // finding node to be deleted
    if (!root) {
        return nullptr;
    }
    if (deletekey < root->key) {
        root->left = remove(root->left, deletekey);
    }
    else if (deletekey > root->key) {
        root->right = remove(root->right, deletekey);
    }
    // this is the node to be deleted
    else {
        // 0 children
        if (!(root->left) && !(root->right)) {
            delete root;
            root = nullptr;
        }
        // 1 child
        else if (!(root->left) || !(root->right)) {
            Node<T>* child;
            if (root->left) {
                child = root->left;
            }
            else {
                child = root->right;
            }
            Node<T>* original = root;
            root = child;
            delete original;
        }
        // 2 children: rotate w/ smaller priority child until it's a leaf
        else {
            if (root->left->priority < root->right->priority) {
                left_rotation(root);
                remove(root->left, deletekey);
            }
            else {
                right_rotation(root);
                remove(root->right, deletekey);
            }
        }
    }
    return root;
}
