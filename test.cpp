
#include <cstdlib>
#include <cstring>
#include <climits>
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

unsigned int getMaxRangeInSet(vector<unsigned int> s)
{
    if (s.size() == 1)
    {
        return s[0];
    }
    else if (s.size() == 2)
    {
        unsigned int m = (s[0] > s[1]) ? s[0] - s[1] : s[1] - s[0];
        return max(m, min(s[0], s[1]));
    }
    else if (s.size() == 3)
    {
        unsigned int x[] = {s[0], s[1], s[2]};
        if (s[1] < s[0]) {
            x[0] = s[1];
            x[1] = s[0];
        }
        if (s[2] < x[1]) {
            x[2] = x[1];
            x[1] = s[2];
            if (x[1] < x[0]) {
                unsigned int temp = x[0];
                x[0] = x[1];
                x[1] = temp;
            }
        }
        unsigned int m = max(x[1] - x[0], x[2] - x[1]);
        return max(m, x[0]);
    }
    else
    {
        unsigned int i, n = s.size(), mn = UINT_MAX, mx = 0, MINVAL = 0, MAXVAL = UINT_MAX;
        for (i = 0; i < n; i++) {
            if (s[i] < mn) {
                mn = s[i];
            }
            if (s[i] > mx) {
                mx = s[i];
            }
        }

        unsigned int * minima_bins = new unsigned int[n - 1];
        unsigned int * maxima_bins = new unsigned int[n - 1];
        for (i = 0; i < n - 1; i++) {
            minima_bins[i] = MAXVAL;
            maxima_bins[i] = MINVAL;
        }

        unsigned int bin_idx;
        double bin_size = (double)(mx - 1) / (double)(n - 1);

        for (i = 0; i < n; i++)
        {
            if (s[i] == mn || s[i] == mx) {
                continue;
            }
            bin_idx = (unsigned int) ((double)(s[i] - mn) / bin_size);
            if (minima_bins[bin_idx] == MAXVAL || s[i] < minima_bins[bin_idx]) {
                minima_bins[bin_idx] = s[i];
            }
            if (maxima_bins[bin_idx] == MINVAL || s[i] > maxima_bins[bin_idx]) {
                maxima_bins[bin_idx] = s[i];
            }
        }

        unsigned int prev = mn, max_gap = 0;
        for (i = 0; i < n - 1; i++)
        {
            if (minima_bins[i] == MAXVAL) {
                continue;
            }
            if ((minima_bins[i] - prev) > max_gap) {
                max_gap = minima_bins[i] - prev;
            }
            prev = maxima_bins[i];
        }

        delete [] minima_bins;
        delete [] maxima_bins;

        if ((mx - prev) > max_gap) {
            max_gap = mx - prev;
        }
        if (mn > max_gap) {
            max_gap = mn;
        }

        return max_gap;
    }
}

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

    cout << endl << endl;

    //testing max range function
    vector<unsigned int> s1 = {33};
    cout << "1: " << getMaxRangeInSet(s1) << endl;
    vector<unsigned int> s2_0 = {15, 18};
    cout << "2.0: " << getMaxRangeInSet(s2_0) << endl;
    vector<unsigned int> s2_1 = {15, 38};
    cout << "2.1: " << getMaxRangeInSet(s2_1) << endl;
    vector<unsigned int> v3_0 = {15, 18, 30};
    vector<unsigned int> v3_1 = {18, 15, 31};
    vector<unsigned int> v3_2 = {2000, 1, 500};
    cout << "3.0: " << getMaxRangeInSet(v3_0) << endl;
    cout << "3.1: " << getMaxRangeInSet(v3_1) << endl;
    cout << "3.2: " << getMaxRangeInSet(v3_2) << endl;
    vector<unsigned int> v4_0 = {15, 18, 30, 10};
    vector<unsigned int> v4_1 = {18, 15, 31, 20};
    vector<unsigned int> v4_2 = {2000, 1, 500, 1599};
    cout << "4.0: " << getMaxRangeInSet(v4_0) << endl;
    cout << "4.1: " << getMaxRangeInSet(v4_1) << endl;
    cout << "4.2: " << getMaxRangeInSet(v4_2) << endl;


    return 0;
}
