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


#endif
