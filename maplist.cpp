#include "hashable.h"
#include "cloneable.h"
#include "map.h"
#include "hashmap.h"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

/**
 * @brief A class which is a wrapper around a vector,
 * but which is also cloneable, and which can therefore
 * be used in our Map interface
 * 
 */
class ValueList: public Cloneable {
    private:
        vector<int> wraplist;
    public:
        ValueList(){};
        ~ValueList(){};
        void push_back(int x) {
            wraplist.push_back(x);
        }
        int get(int i) {
            return wraplist[i];
        }
        size_t size() {
            return wraplist.size();
        }
        Cloneable* clone() {
            ValueList* ret = new ValueList();
            for (size_t i = 0; i < size(); i++) {
                ret->push_back(wraplist[i]);
            }
            return ret;
        }
        void print() {
            for (size_t i = 0; i < size(); i++) {
                cout << wraplist[i];
                if (i < size()-1) {
                    cout << ", ";
                }
            }
            cout << "\n";
        }
};

int main(int argc, char** argv) {
    Map* m = new HashMap(100);
    
    // Here's how to add the lists to the map, using
    // a hashable key (a HashableString in this case)
    HashableString s1(string("list1"));
    ValueList list1;
    list1.push_back(0);
    list1.push_back(10);
    m->put(&s1, &list1);

    HashableString s2(string("list2"));
    ValueList list2;
    m->put(&s2, &list2);

    // Here's how to get the list from the map and add stuff to it
    HashableString query(string("list1"));
    ValueList* list = (ValueList*)m->get(&query);
    list->push_back(20);
    list->push_back(30);
    list->print(); // This should output 0, 10, 20, 30

    delete m;
    return 0;
}

