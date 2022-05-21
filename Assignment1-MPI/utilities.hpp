// ==============================================================================
// 
//  Write your code in "gensupermer_mpi.cpp" not this file. You may modify
//  this file for debugging but anything wrote to this file will not be counted.
//
//  You should NOT use any variable or funtion from this file in your section.
// 
// ==============================================================================
// 7 MAR (update2.1): update GenerateSupermer_seq, 28 FEB: update line #27

#pragma once
#define _in_
#define _out_

#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
using namespace std;

static char* getline(FILE* fp, char *buffer, /**/int buffer_size, _out_ string &line) {
    line.clear();
    char* res_flag = fgets(buffer, buffer_size-1, fp);
    char* t = res_flag;
    while (buffer[strlen(buffer)-1] != '\r' && buffer[strlen(buffer)-1] != '\n' && t && !feof(fp)) {
        line += buffer;
        t = fgets(buffer, buffer_size-1, fp);
    } 
    for (int i = strlen(buffer)-1; i>=0; i--) {
        if (buffer[i] == '\r' || buffer[i] == '\n') buffer[i] = 0;
        else break;
    }
    line += buffer;
    return res_flag;
}

static void LoadReadsFromFile(const char* filename, _out_ vector<string> &reads) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        cerr << "Error when open " << filename << "." << endl;
        exit(1);
    }
    const int BUF_SIZE = 1048576 * 10;
    string line;
    char *read_buffer = new char[BUF_SIZE];
    while(getline(fp, read_buffer, BUF_SIZE, line)) {
        if (line.length()>22) reads.push_back(line);
    }
    fclose(fp);
    return;
}

/**
* @brief CSR Format
* @par Sample:
* @code
*   //sample code here
* @endcode
* @note Limitations:
* The length of one single data item cannot exceed INT_MAX;
* The total number of items cannot exceed INT_MAX;
* The total length of the whole data cannot exceed ULL_MAX;
*/
typedef /**/int T_CSR_capacity;
typedef int T_CSR_count;
template<typename T_data, class T_attr=bool> class CSR
{
protected:
    T_data *_data;
    bool _with_attr;
    T_attr *_attr;
    T_CSR_count _item_count = 0, _item_capacity = 4; // for saving attributes
    T_CSR_capacity *_offs;
    
    T_CSR_capacity _size = 0;     // unit: counts not bytes, for saving raw data
    T_CSR_capacity _capacity = 0; // unit: counts not bytes
    
    void _Realloc(T_CSR_capacity new_capacity) {
        // new_capacity: capacity in counts
        _data = (T_data*)realloc(_data, new_capacity * sizeof(T_data));
    }
    
public:
    float size_incre_ratio = 1.5;
    
    CSR(bool with_attr = false, T_CSR_capacity initial_capacity = 8) {
        _with_attr = with_attr;
        if (with_attr) _attr = new T_attr[_item_capacity]();//
        _offs = new T_CSR_capacity[_item_capacity+1]();//
        _offs[0] = 0;
        _capacity = initial_capacity;
        _data = new T_data[_capacity]();//
    }
    CSR(const CSR<T_data, T_attr>& f) {
        _item_capacity = _item_count = f._item_count;
        if (f.has_attr()) {
            _attr = new T_attr[_item_capacity]();//
            _with_attr = true;
            memcpy(_attr, f._attr, sizeof(T_attr) * _item_count);
        }
        _offs = new T_CSR_capacity[_item_capacity+1]();//?
        memcpy(_offs, f._offs, sizeof(T_CSR_capacity)*(_item_count+1));
        _size = _capacity = f._size;
        _data = new T_data[_size]();//
        memcpy(_data, f._data, sizeof(T_data)*_size);
    }
    ~CSR(void) {
        delete [] _offs;//
        delete [] _data;//
        if (_with_attr) delete [] _attr;//
    }
    
    // append with no attr
    void append(const T_data *new_data, T_CSR_capacity new_data_size) {
        // the unit of new_data_size is count not byte
        if (new_data_size + _size >= _capacity) {
            T_CSR_capacity new_capacity = (T_CSR_capacity)(double(_capacity)*size_incre_ratio);
            _capacity = new_data_size + _size > new_capacity ? new_data_size + _size + 1 : new_capacity;
            _Realloc(_capacity);
        }
        if (_item_count == _item_capacity-1) {
            _offs = (T_CSR_capacity *)realloc(_offs, sizeof(T_CSR_capacity) * (2*_item_capacity+1));
            if (_with_attr) _attr = (T_attr *)realloc(_attr, sizeof(T_attr) * (2*_item_capacity));
            _item_capacity = 2*_item_capacity;
        }
        memcpy(_data + _offs[_item_count], new_data, new_data_size * sizeof(T_data));
        _size += new_data_size;
        _offs[_item_count+1] = _size;
        _item_count += 1;
    }

    // append with attr
    void append(const T_data *new_data, T_CSR_capacity new_data_size, T_attr attr) {
        // the unit of new_data_size is count not byte
        if (new_data_size + _size >= _capacity) {
            T_CSR_capacity new_capacity = (T_CSR_capacity)(double(_capacity)*size_incre_ratio);
            _capacity = new_data_size + _size > new_capacity ? new_data_size + _size + 1 : new_capacity;
            _Realloc(_capacity);
        }
        if (_item_count == _item_capacity-1) {
            _offs = (T_CSR_capacity *)realloc(_offs, sizeof(T_CSR_capacity) * (2*_item_capacity+1));
            if (_with_attr) _attr = (T_attr *)realloc(_attr, sizeof(T_attr) * (2*_item_capacity));
            _item_capacity = 2*_item_capacity;
        }
        memcpy(_data + _offs[_item_count], new_data, new_data_size * sizeof(T_data));
        _size += new_data_size;
        if (_with_attr) _attr[_item_count] = attr;
        _offs[_item_count+1] = _size;
        _item_count += 1;
    }
    void append(CSR<T_data, T_attr> new_data) {
        T_CSR_capacity new_data_size = new_data.size();
        if (new_data_size + _size >= _capacity) {
            T_CSR_capacity new_capacity = (T_CSR_capacity)(double(_capacity)*size_incre_ratio);
            _capacity = new_data_size + _size > new_capacity ? new_data_size + _size + 1 : new_capacity;
            _Realloc(_capacity);
        }
        if (_item_count + new_data.items() >= _item_capacity) {
            // should multiply 2 below to avoid float size.
            int new_size = _item_count+new_data.items() > 2*_item_capacity ? _item_count+new_data.items() : 2*_item_capacity;
            _offs = (T_CSR_capacity *)realloc(_offs, sizeof(T_CSR_capacity) * (new_size+1));
            if (_with_attr) _attr = (T_attr *)realloc(_attr, sizeof(T_attr) * (new_size));
            _item_capacity = new_size;
        }
        memcpy(_data + _offs[_item_count], new_data.get_raw_data(), new_data_size * sizeof(T_data));
        memcpy(&_offs[_item_count+1], new_data.get_raw_offs()+1, new_data.items() * sizeof(T_CSR_capacity));
        for(T_CSR_count i=_item_count+1; i<_item_count+new_data.items()+1; i++)
            _offs[i] += _offs[_item_count]; // add the offset value for new elements
        if (_with_attr && new_data.has_attr()) 
            memcpy(&_attr[_item_count], new_data.get_raw_attr(), new_data.items() * sizeof(T_attr));
        _size += new_data_size;
        _item_count += new_data.items();
    }
    T_CSR_capacity capacity() {
        return _capacity;
    }
    T_CSR_capacity size() {
        return _size;
    }
    T_CSR_count items() {
        return _item_count;
    }
    short dtype() {
        return sizeof(T_data);
    }
    void debug_info() {
        cout << "ITEMS\t" << this -> items() << "/" << _item_capacity << endl;
        for (int i=0; i<_item_count; i++) {
            cout << "offs="<< _offs[i+1];
            if (has_attr()) cout << " attr=" << _attr[i];
            cout << endl;
        }
        cout << "SIZE\t" << this -> size() << "/" << this -> capacity() << endl;
        cout << "DTYPE\t" <<this -> dtype() << endl;
        for (T_CSR_capacity i=0; i<_size; i++) cout<<_data[i]<<" ";
        cout << endl;
    }
    T_data* get_raw_data() {
        return _data;
    }
    T_CSR_capacity* get_raw_offs() {
        return _offs;
    }
    T_data* fetch_raw_data() {
        T_data* res = new T_data[_size];
        memcpy(res, _data, _size*sizeof(T_data));
        return res;
    }
    T_CSR_capacity* fetch_raw_offs() {
        T_CSR_capacity* res = new T_CSR_capacity[_item_count+1];
        memcpy(res, _offs, (_item_count+1)*sizeof(T_CSR_capacity));
        return res;
    }
    bool has_attr() {
        return _with_attr;
    }
    T_attr* get_raw_attr() {
        return _attr;
    }
    int get_item(int idx, _out_ T_data *data_buffer) {
        if (idx > this->size()) return -1;
        memcpy(data_buffer, &(this->_data[this->_offs[idx]]), sizeof(T_data) * (this->_offs[idx+1] - this->_offs[idx]));
        return 0;
    }
    T_attr get_attr(int idx) {
        if (idx > this->size()) {
            cerr << "[ERROR] CSR.get_attr(idx) got wrong idx " << idx << endl;
            exit(1);
        }
        return _attr[idx];
    }
};

static void Vector2CSR(vector<string> &reads, int &num_of_reads, char* &reads_CSR, /**/int* &reads_CSR_offs) {
    CSR<char> t;
    for(string read: reads) {
        t.append(read.c_str(), read.length());
    }
    num_of_reads = reads.size();
    reads_CSR = t.fetch_raw_data();
    reads_CSR_offs = t.fetch_raw_offs();
}

void GenerateSupermer_seq(vector<string> reads, int k, int p, _out_ vector<string> &supermers) { // UPDATE 2.1
    for (string read: reads) {
        int read_len = read.length();
        string prev_minimizer, minimizer, new_minimizer;
        int i, j;
        char base;
        int skm_begin_pos, skm_end_pos, mm_begin_pos;

        // Generate the first k-mer's minimizer:
        skm_begin_pos = 0;
        skm_end_pos = k;
        mm_begin_pos = 0;
        minimizer = new_minimizer = read.substr(0, p);
        for (i=p; i<k; i++) {
            new_minimizer = new_minimizer.substr(1, p-1) + read[i]; // UPDATE1
            if (new_minimizer <= minimizer) minimizer = new_minimizer, mm_begin_pos = i-p+1;
        }

        // Continue generating minimizers:
        for (i=1; i<read_len-k+1; i++) { // i: the beginning position of the current k-mer
            if (i > mm_begin_pos) {
                // new minimizer required
                prev_minimizer = minimizer;
                minimizer = new_minimizer = read.substr(i, p);
                for (j=i+p; j<i+k; j++) {
                    new_minimizer = new_minimizer.substr(1, p-1) + read[j]; // UPDATE1
                    if (new_minimizer <= minimizer) minimizer = new_minimizer, mm_begin_pos = j-p+1;
                }
                // if the new minimizer equals to the previous one, we can continue
                if (minimizer != prev_minimizer) {
                    skm_end_pos = i-1+k;
                    supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos)); // save the supermer
                    skm_begin_pos = i;
                }
            }
            else {
                new_minimizer = read.substr(i+k-p, p); // UPDATE1
                if (new_minimizer < minimizer) { // save the supermer
                    skm_end_pos = i-1+k;
                    supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos));
                    skm_begin_pos = i;
                    minimizer = new_minimizer, mm_begin_pos = i+k-1-p+1;
                }
                if (new_minimizer == minimizer) mm_begin_pos = i+k-1-p+1; // UPDATE1
            }
        } // UPDATE 2.1
        skm_end_pos = read_len;
        supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos));
    }
}

static void CorrectnessChecking(vector<string> &reads, int k, int p, vector<string> &supermers) {
    vector<string> supermers_gt;
    cout << "Now using sequential algorithm to check correctness (may take some time) ..." << endl;
    auto start_time = chrono::high_resolution_clock::now();
    GenerateSupermer_seq(reads, k, p, supermers_gt);
    auto end_time = chrono::high_resolution_clock::now(); // time measurement ends
    auto duration_sec = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count()/1000.0;
    cout << "Sequential-version supermer generation ends in " << duration_sec << " sec." << endl << endl; // 173.906s
    sort(supermers.begin(), supermers.end());
    sort(supermers_gt.begin(), supermers_gt.end());
    if (supermers == supermers_gt) {
        cout << "Test passed!" << endl << endl;
    }
    else {
        cout << "Test failed, below is the first 10 incorrect supermers (yours/ground_truth):" << endl;
        int error_count = 0;
        vector<string>::iterator i, j;
        for (i = supermers.begin(), j = supermers_gt.begin(); 
                i != supermers.end() && j != supermers_gt.end() && error_count < 10; ) {
            if (*i != *j) {
                error_count ++;
                if (*i<*j) cout<<"WRONG" << error_count << ":\t" << *i << "\t" << *j << endl;
                else cout<<"MISS " << error_count << ":\t -------- " << "\t" << *j << endl;
                if (*i < *j) i++;
                else j++;
            }
            else i++, j++;
        }
        cout << "Test failed!" << endl << endl;
    }
}

void ArgParser(int &argc, char **&argvs, int &k, int &p, string &read_file_path, int &correctness_checking, string &output_path) {
    if (argc < 5) {
        cerr << "Error! usage:" << endl << "\"./gs [k] [p] [data file path] [correctness_checking=0/1] (optional) [supermers_output_path]\"" << endl;
        exit(1);
    }
    k = atoi(argvs[1]);
    p = atoi(argvs[2]);
    read_file_path = argvs[3];
    if (argc >= 5) correctness_checking = atoi(argvs[4]);
    if (argc >= 6) {
        output_path = string(argvs[5]);
        if (*(output_path.end()-1) != '/') output_path += "/my_gs_output.txt";
        else output_path += "my_gs_output.txt";
    }
}

void SaveSupermers(string &output_path, vector<string> &supermers) {
    FILE *fp = fopen(output_path.c_str(), "w");
    if (fp == NULL) {
        cerr << "ERROR: unable to write " << output_path << ", please check your permission." << endl;
        exit(1);
    }
    for (string supermer: supermers) fprintf(fp, "%s\n", supermer.c_str());
    fclose(fp);
    cout << supermers.size() << " supermers have been saved to " << output_path << "." << endl;
}

// void GenerateSupermer(vector<string> reads, int k, int p, _out_ vector<string> &supermers);

