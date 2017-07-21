
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <queue>
#include <sdsl/util.hpp>
#include <sdsl/suffix_trees.hpp>

typedef sdsl::cst_sct3<> cst_t;
typedef cst_t::node_type cst_node_t;
typedef cst_t::size_type cst_size_t;
typedef cst_t::char_type cst_char_t;

using namespace std;
using namespace sdsl;

int main()
{
    string p = "ACACACACCAACAAC";
    //           01234567890123456
    string sp = "ACACA#CACCA#ACAAC";
    int STpIdx2BVIdx[] = {0,1,2,3,4,-1,5,6,7,8,9,-1,10,11,12,13,14,-1};
    cst_t STp;
    construct_im(STp, sp.c_str(), sizeof(char));

    //print idx, sn, suffix, lcp of csa

    cout << "idx\tsn\tlcp\tsuffix" << endl;
    unsigned int i;
    for (i = 0; i < sp.length() + 1; i++)
    {
        cout << i << "\t" \
             << STp.csa[i] << "\t" \
             << STp.lcp[i] << "\t" \
             << sp.substr(STp.csa[i]) << endl;
    }
    cout << endl;

    //print nodes in tree

    cout << "nid\tlb\trb\tsn\tleaf\ttext" << endl;
    queue<cst_node_t> q;
    q.push(STp.root());
    cst_node_t currNode;
    while (!q.empty())
    {
        currNode = q.front();
        q.pop();

        cout << STp.id(currNode) << "\t" \
             << STp.lb(currNode) << "\t" \
             << STp.rb(currNode) << "\t";
        cout << ((STp.is_leaf(currNode)) ? (int) STp.sn(currNode) : -1) << "\t";
        cout << ((STp.is_leaf(currNode)) ? "Y" : " ") << "\t";
        if (currNode != STp.root()) {
            cout << extract(STp, currNode);
        }
        cout << endl;

        for (const auto & child : STp.children(currNode))
        {
            q.push(child);
        }
    }
    cout << endl << endl;

    //identify which nodes of CA we should follow: 13, 24, 1, 6 but not 23, 9, 3.

    cout << "Traversing down CA..." << endl;
    while (!q.empty()) {
        q.pop();
    }
    q.push(STp.inv_id(25)); //CA
    unsigned int level = 2;
    unsigned int dc = 1; //decrement counter
    unsigned int ic = 0; //increment counter
    while (!q.empty())
    {
        currNode = q.front();
        q.pop();
        for (const auto & child : STp.children(currNode))
        {
            unsigned int id = STp.id(child);
            if (!STp.is_leaf(child))
            {
                cout << "Visited node " << id << " on level " << (level + 1) << endl;
                unsigned int lb = STp.lb(child);
                unsigned int sn = STp.csa[lb];
                char nextChar = sp[sn + level];
                if (nextChar != '#') {
                    q.push(child);
                }
            }
            else
            {
                cout << "Visited child " << id << " on level " << (level + 1) << endl;
            }
            ic++;
        }
        dc--;
        if (dc == 0) {
            level++;
            dc = ic;
            ic = 0;
        }
    }
    cout << endl << endl;

    //identify which nodes of A we should follow: 14, 21, 20, 7, 2, 12, 0 but not 10, 4

    cout << "Traversing down A..." << endl;
    while (!q.empty()) {
        q.pop();
    }
    q.push(STp.inv_id(22)); //A
    level = 1;
    dc = 1; //decrement counter
    ic = 0; //increment counter
    while (!q.empty())
    {
        currNode = q.front();
        q.pop();
        for (const auto & child : STp.children(currNode))
        {
            unsigned int id = STp.id(child);
            if (!STp.is_leaf(child))
            {
                cout << "Visited node " << id << " on level " << (level + 1) << endl;
                unsigned int lb = STp.lb(child);
                unsigned int sn = STp.csa[lb];
                char nextChar = sp[sn + level];
                if (nextChar != '#') {
                    q.push(child);
                }
            }
            else
            {
                cout << "Visited child " << id << " on level " << (level + 1) << endl;
            }
            ic++;
        }
        dc--;
        if (dc == 0) {
            level++;
            dc = ic;
            ic = 0;
        }
    }

    return 0;
}
