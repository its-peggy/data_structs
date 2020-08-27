/*
Sep 2019
Implement Mod class
*/

#include "mod.h"

Mod::Mod(long t) {
    long res = t % modulus;
    if (res < 0) { x = (res + modulus) % modulus; }
    else { x = res; }
}

Mod::Mod(const Mod& m) : x(m.x) {}

Mod& Mod::operator=(const Mod& m) {
    x = m.x;
    return *this;
}

Mod& Mod::operator+=(const Mod& m) {
    long a = x;
    long b = m.x;
    if (a < 0) { a += modulus; }
    if (b < 0) { b += modulus; }
    long res = (a - modulus) + b;
    if (res < 0) { res += modulus; }
    x = res;
    return *this;
}

Mod& Mod::operator-=(const Mod& m) {
    long a = x;
    long b = m.x;
    if (a < 0) { a += modulus; }
    if (b < 0) { b += modulus; }
    long res = a - b;
    if (res < 0) { res += modulus; }
    x = res;
    return *this;
}

Mod& Mod::operator*=(const Mod& m) {
    Mod cur = m;
    Mod res = Mod(0);
    while (x > 0) {
        int bit = x & 1;
        if (bit == 1) {
            res += cur;
        }
        cur += cur;
        x = (x >> 1);
    }
    x = res.x;
    return *this;
}

Mod& Mod::operator/=(const Mod& m) {
    Mod res = Mod(x);
    res *= inv(m.x);
    x = res.x;
    return *this;
}

Mod Mod::operator-() const {
    Mod y = Mod(x);
    Mod z = Mod(modulus-y.x);
    return z;
}

Mod Mod::pwr(long e) const {
    if (e < 0) {
        return inv(x).pwr(-1 * e);
    }
    else {
        Mod cur = x;
        Mod res = Mod(1);
        while (e > 0) {
            int bit = e & 1;
            if (bit == 1) {
                res *= cur;
            }
            cur *= cur;
            e = (e >> 1);
        }
        return res;
    }
}

long Mod::val() const {
    return x;
}

void Mod::set_modulus(long m) {
    modulus = m;
}

long x;

long modulus = 17;

Mod Mod::inv(long r0) {
    r0 = (r0 % modulus);
    if (r0 < 0) {
        r0 = (r0 + modulus) % modulus;
    }

    long j = 0;
    long k = 1;
    long q = modulus;
    long old_j = 1;
    long old_k = 0;
    long old_q = r0;

    while (q > 0) {
        long div = old_q / q;
        long temp;
        temp = q;
        q = old_q - div * q;
        old_q = temp;
        temp = j;
        j = old_j - div * j;
        old_j = temp;
        temp = k;
        k = old_k - div * k;
        old_k = temp;
    }
    if (old_q != 1) {
        cerr << "param and modulus not relatively prime...\n";
        exit(-1);
    }
    Mod res = Mod(old_j);
    return res;
}

//////////////////////////////////////////////////////////////////////

Mod operator+(const Mod& a, const Mod& b) {
    Mod res = Mod(a.val());
    res += b;
    return res;
}

Mod operator+(long t, const Mod& m) {
    Mod res = Mod(t);
    res += m;
    return res;
}

Mod operator-(const Mod& a, const Mod& b) {
    Mod res = Mod(a.val());
    res -= b;
    return res;
}

Mod operator-(long t, const Mod& m) {
    Mod res = Mod(t);
    res -= m;
    return res;
}

Mod operator*(const Mod& a, const Mod& b) {
    Mod res = Mod(a.val());
    res *= b;
    return res;
}

Mod operator*(long t, const Mod& m) {
    Mod res = Mod(t);
    res *= m;
    return res;
}

Mod operator/(const Mod& a, const Mod& b) {
    Mod res = Mod(a.val());
    res /= b;
    return res;
}

Mod operator/(long t, const Mod& m) {
    Mod res = Mod(t);
    res /= m;
    return res;
}

bool operator==(const Mod& a, const Mod& b) {
    return (a.val() == b.val());
}

bool operator==(long t, const Mod& m) {
    Mod m2 = Mod(t);
    return (m2.val() == m.val());
}

bool operator!=(const Mod& a, const Mod& b) {
    return (a.val() != b.val());
}

bool operator!=(long t, const Mod& m) {
    Mod m2 = Mod(t);
    return (m2.val() != m.val());
}

istream& operator>>(istream& is, Mod& m) {
    long val;
    is >> val;
    m = Mod(val);
    return is;
}

ostream& operator<<(ostream& os, const Mod& m) {
    os << m.val();
    return os;
}

int main() {
    cout << "\n";
    Mod::set_modulus(101);
    long mod = Mod::get_modulus();
    cout << "modulus: " << mod << "\n";
    Mod one(1);
    Mod m(2);
    cout << one/m << "\n";
    //cout << one/cur << "\n";
}
