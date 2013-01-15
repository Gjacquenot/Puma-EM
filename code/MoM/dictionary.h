#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <vector>
#include <algorithm>

// we need the following Dictionary class in various places
template <typename K, typename V>
class Dictionary {
    K key;
    V val;
  public:
    // constructors
    Dictionary(const K k, const V v) {key = k; val = v;};
    const V getVal(void) const {return val;};
    const K getKey(void) const {return key;};
    Dictionary(const Dictionary & dictionaryToCopy) {val = dictionaryToCopy.getVal(); key = dictionaryToCopy.getKey();}; // copy constructor
    ~Dictionary(){};
    // overloaded operators for sorting
    bool operator== (const Dictionary & right) const {if ( getKey() == right.getKey() ) return 1; else return 0;};
    bool operator< (const Dictionary & right) const {if ( getKey() < right.getKey() ) return 1; else return 0;};
};

template <typename K, typename V1, typename V2>
class Dictionary2 {
    K key;
    V1 val1;
    V2 val2;
  public:
    // constructors
    Dictionary2(const K k, const V1 v1, const V2 v2) {key = k; val1 = v1; val2 = v2;};
    const V1 getVal1(void) const {return val1;};
    const V2 getVal2(void) const {return val2;};
    const K getKey(void) const {return key;};
    Dictionary2(const Dictionary2 & dictionary2ToCopy) {val1 = dictionary2ToCopy.getVal1(); val2 = dictionary2ToCopy.getVal2(); key = dictionary2ToCopy.getKey();}; // copy constructor
    ~Dictionary2(){};
    // overloaded operators for sorting
    bool operator== (const Dictionary2 & right) const {if ( getKey() == right.getKey() ) return 1; else return 0;};
    bool operator< (const Dictionary2 & right) const {if ( getKey() < right.getKey() ) return 1; else return 0;};
};

#endif
