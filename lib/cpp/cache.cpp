#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

class LFUCache {
public:
    LFUCache(int capacity) {
        _capacity = capacity;
        _min_freq = 0;
    }

    int get(int key) {
        if (_cache.count(key) == 0) {
            return -1;
        }

        // increase frequency of key
        int freq = _cache[key].second++;
        _freq_list[freq].erase(_key_list[key]);
        _freq_list[freq + 1].push_back(key);
        _key_list[key] = --_freq_list[freq + 1].end();

        // update min_freq
        if (_freq_list[_min_freq].size() == 0) {
            _min_freq++;
        }

        return _cache[key].first;
    }

    void put(int key, int value) {
        if (_capacity <= 0) {
            return;
        }

        if (get(key) != -1) {
            _cache[key].first = value;
            return;
        }

        if (_cache.size() >= _capacity) {
            // remove least frequently used key
            int least_freq_key = _freq_list[_min_freq].front();
            _freq_list[_min_freq].pop_front();
            _cache.erase(least_freq_key);
            _key_list.erase(least_freq_key);
        }

        // add new key
        _cache[key] = make_pair(value, 1);
        _freq_list[1].push_back(key);
        _key_list[key] = --_freq_list[1].end();
        _min_freq = 1;
    }

private:
    int _capacity;
    int _min_freq;
    map<int, pair<int, int>> _cache; // key -> (value, frequency)
    map<int, list<int>> _freq_list; // frequency -> key list
    map<int, list<int>::iterator> _key_list; // key -> iterator in freq_list
};

class LRUCache {
public:
    LRUCache(int capacity) {
        _capacity = capacity;
    }

    int get(int key) {
        auto it = _cache.find(key);
        if (it == _cache.end()) {
            return -1;
        }

        // move the key to the front of the list
        _key_list.splice(_key_list.begin(), _key_list, it->second);
        return it->second->second;
    }

    void put(int key, int value) {
        auto it = _cache.find(key);
        if (it != _cache.end()) {
            // update the value
            it->second->second = value;
            // move the key to the front of the list
            _key_list.splice(_key_list.begin(), _key_list, it->second);
            return;
        }

        if (_cache.size() == _capacity) {
            // evict the least recently used key
            int least_used_key = _key_list.back().first;
            _cache.erase(least_used_key);
            _key_list.pop_back();
        }

        // insert the new key-value pair
        _key_list.emplace_front(key, value);
        _cache[key] = _key_list.begin();
    }

private:
    int _capacity;
    list<pair<int, int>> _key_list;
    unordered_map<int, list<pair<int, int>>::iterator> _cache;
};

int main() {
  return 0;
}
