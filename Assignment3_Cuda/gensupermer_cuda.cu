//[NAME]: QIU Yaowen
//[STUDENT ID]:   20784389
//[EMAIL]:    yqiuau@connect.ust.hk

#include <iostream>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <assert.h>
#include <thread>
#include <future>
#include <string>
#include "gensupermer.hpp"

#include "utilities.hpp"
using namespace std;

__device__ __constant__ static const unsigned char d_basemap[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, // 0..7
    255, 255, 255, 255, 255, 255, 255, 255, // 8..15
    255, 255, 255, 255, 255, 255, 255, 255, // 16..23
    255, 255, 255, 255, 255, 255, 255, 255, // 24..31
    255, 255, 255, 255, 255, 255, 255, 255, // 32..39
    255, 255, 255, 255, 255, 255, 255, 255, // 40..47
    255, 255, 255, 255, 255, 255, 255, 255, // 48..55
    255, 255, 255, 255, 255, 255, 255, 255, // 56..63
    255, 0, 255, 1, 255, 255, 255, 2, // 64..71
    255, 255, 255, 255, 255, 255, 255, 255, // 72..79
    255, 255, 255, 255, 3, 255, 255, 255, // 80..87
    255, 255, 255, 255, 255, 255, 255, 255, // 88..95
    255, 0, 255, 1, 255, 255, 255, 2, // 96..103
    255, 255, 255, 255, 255, 255, 255, 255, // 104..111
    255, 255, 255, 255, 3, 255, 255, 255, // 112..119
    255, 255, 255, 255, 255, 255, 255, 255, // 120..127
    255, 255, 255, 255, 255, 255, 255, 255, // 128..135
    255, 255, 255, 255, 255, 255, 255, 255, // 136..143
    255, 255, 255, 255, 255, 255, 255, 255, // 144..151
    255, 255, 255, 255, 255, 255, 255, 255, // 152..159
    255, 255, 255, 255, 255, 255, 255, 255, // 160..167
    255, 255, 255, 255, 255, 255, 255, 255, // 168..175
    255, 255, 255, 255, 255, 255, 255, 255, // 176..183
    255, 255, 255, 255, 255, 255, 255, 255, // 184..191
    255, 255, 255, 255, 255, 255, 255, 255, // 192..199
    255, 255, 255, 255, 255, 255, 255, 255, // 200..207
    255, 255, 255, 255, 255, 255, 255, 255, // 208..215
    255, 255, 255, 255, 255, 255, 255, 255, // 216..223
    255, 255, 255, 255, 255, 255, 255, 255, // 224..231
    255, 255, 255, 255, 255, 255, 255, 255, // 232..239
    255, 255, 255, 255, 255, 255, 255, 255, // 240..247
    255, 255, 255, 255, 255, 255, 255, 255  // 248..255
};

typedef struct {
    _in_ T_read_count cur_batch_size;    //number of reads
    
    // Raw reads
    _in_ _out_ char *reads;              //read_CSR 
    _in_ T_CSR_capacity *reads_offs;     //read_CSR_offset 
    _in_ _out_ T_read_len *read_len;      //length of each read
    
    // Minimizers
    _out_ _tmp_ T_minimizer *minimizers;  //the minimizer array
    _out_ T_read_len *supermer_offs;     //the supermer offset array
} T_GPU_data;

/*
 * [INPUT]  data.reads in [(Read#0), (Read#1)...]
 * [OUTPUT] data.minimizers in [(Read#0)[mm1, mm?, mm?, ...], (Read#1)...]
 */
__global__ void GPU_GenMinimizer(_in_ _out_ T_GPU_data data, int K_kmer, int P_minimizer) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x; //index of current kernel
    const int nthread = blockDim.x * gridDim.x;

    int index = tid;
    int batch_size = data.cur_batch_size;
    // if (index == 0){
    //     printf("%d %d\n",data.reads_offs[index],data.read_len[index]);
    //     printf("%d %d\n",data.reads_offs[4095],data.read_len[4095]);
    //     printf("%d %d\n",data.reads_offs[4096],data.read_len[4096]);
    //     printf("%d %d\n",data.reads_offs[4400],data.read_len[4400]);
    // }
    // for(int index = tid; index < batch_size; index += nthread) {
    while (true){
        // printf("Index :%d,offset:%d \n",index,data.reads_offs[index]);

        if (index >= batch_size){
            return;
        }
        
        // if (index != tid){
        //     printf("Index : %d, Len : %d \n",index,data.read_len[index]);
        // }

        // if (data.read_len[index] == 0){
        //     return;
        // }

        int n = data.read_len[index]; //length of current read

        unsigned int startIdx = data.reads_offs[index];
        unsigned int insert_index = startIdx; //indicate where current minimizer should insert into data.minimizers
        int counter = 0; //counter how many insertion has made

        unsigned int minimizer,new_minimizer,prev_minimizer;
        minimizer = new_minimizer = prev_minimizer = 0;
        unsigned int ptr_new_minimizer = 0 + startIdx; //pointer indicate the start position of new_minimizer
        unsigned int tmp_minimizer = 0; //tmp varaible to get substr representation

        int i,j,m;
        i = j = m = 0;

        unsigned int skm_begin_pos, skm_end_pos, mm_begin_pos;

        // **** Generate the first k-mer's minimizer Starts
        skm_begin_pos = 0 + startIdx;
        skm_end_pos = K_kmer + startIdx;
        mm_begin_pos = 0 + startIdx;

        for(i = 0 + startIdx;i<P_minimizer + startIdx;i++){
            minimizer = (minimizer << 2) | d_basemap[data.reads[i]]; // minimizer assignment
        }


        new_minimizer = minimizer; // same as minimizer

        for (i=P_minimizer + startIdx; i<K_kmer + startIdx; i++) {

            // recompute new_minimizer
            new_minimizer = 0;
            ptr_new_minimizer += 1;
            for (m = ptr_new_minimizer; m < ptr_new_minimizer + P_minimizer - 1; m++){
                new_minimizer = (new_minimizer << 2) | d_basemap[data.reads[m]]; //simulate sliding windows
            }
            new_minimizer = (new_minimizer << 2) | d_basemap[data.reads[i]];

            if (new_minimizer <= minimizer) {
                minimizer = new_minimizer;
                mm_begin_pos = i-P_minimizer+1;
            }
        }

        data.minimizers[insert_index] = minimizer;
        insert_index ++;
        counter ++;
        // **** Generate the first k-mer's minimizer Ends

        // ** Continue generating minimizers Starts
        for (i = 1 + startIdx; i < n- K_kmer + 1 + startIdx; i++) { // i: the beginning position of the current k-mer
            if (i > mm_begin_pos) {

                // new minimizer required
                prev_minimizer = minimizer;
                
                tmp_minimizer = 0;
                for (m = i; m < i + P_minimizer; m++){
                    tmp_minimizer = (tmp_minimizer << 2) | d_basemap[data.reads[m]];
                }
                
                minimizer = new_minimizer = tmp_minimizer;
                ptr_new_minimizer = i; //update pointer of new_minimizer
                tmp_minimizer = 0;

                for (j=i+P_minimizer; j<i+K_kmer; j++) {

                    // recompute new_minimizer
                    new_minimizer = 0;
                    ptr_new_minimizer += 1;
                    for (m = ptr_new_minimizer; m < ptr_new_minimizer + P_minimizer-1;m++){
                        new_minimizer = (new_minimizer << 2) | d_basemap[data.reads[m]]; //simulate sliding windows
                    }
                    new_minimizer = (new_minimizer << 2) | d_basemap[data.reads[j]];
                    
                    if (new_minimizer <= minimizer){
                        minimizer = new_minimizer;
                        mm_begin_pos = j-P_minimizer+1;
                    }
                }

                // if the new minimizer equals to the previous one, we can continue
                if (minimizer != prev_minimizer) {
                    skm_end_pos = i-1+K_kmer;
                    skm_begin_pos = i;
                }
            }
            else {

                tmp_minimizer = 0;
                for(m = i + K_kmer - P_minimizer; m < i + K_kmer; m++){
                    tmp_minimizer = (tmp_minimizer << 2) | d_basemap[data.reads[m]];
                }
                new_minimizer = tmp_minimizer; // UPDATE1
                ptr_new_minimizer = i + K_kmer - P_minimizer;
                tmp_minimizer = 0;

                if (new_minimizer < minimizer) { // save the supermer
                    skm_end_pos = i-1+K_kmer;
                    skm_begin_pos = i;
                    minimizer = new_minimizer;
                    mm_begin_pos = i+K_kmer-1-P_minimizer+1;
                }
                if (new_minimizer == minimizer){
                    mm_begin_pos = i+K_kmer-1-P_minimizer+1; // UPDATE1
                }
            }

            data.minimizers[insert_index] = minimizer;
            insert_index ++;
            counter ++;

        }
        // UPDATE 2.1
        // ** Continue generating minimizers Ends


        while (counter < n ){ //fill in null value
            data.minimizers[insert_index] = 88888888; // insert nothing to the output array
            insert_index++;
            counter++;
        }

        // if (index == 100){
        //     for(int i = data.reads_offs[index]; i < data.reads_offs[index] + n; i++){
                
        //         printf("Index: %d, Minimzer: %d ",i - data.reads_offs[index],data.minimizers[i]);
                
        //     }
        // }

        
        index += nthread;
        
        // return;
    }
    
}


/* [INPUT]  data.minimizers in [[mm1, mm1, mm2, mm3, ...], ...]
 * [OUTPUT] data.supermer_offs in [[0, 2, 3, ...], ...]
 */
__global__ void GPU_GenSKM(_in_ _out_ T_GPU_data data, int K_kmer, int P_minimizer) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x; //index of current kernel
    const int nthread = blockDim.x * gridDim.x;

    int index = tid;
    int batch_size = data.cur_batch_size;
    // for(int index = tid; index < batch_size; index += nthread) {
    while (true){

        if (index >= batch_size){
            return;
        }

        // if (data.read_len[index] == 0){
        //     return;
        // }
        
        int n = data.read_len[index]; //length of current read
        unsigned int startIdx = data.reads_offs[index];
        unsigned long insert_index = startIdx; //indicate where current supermer_offs should insert into data.supermer_offs

        // if (index == 100){
        //     for (int j = startIdx; j < startIdx + 10;j++){
        //         printf("Index : %d, Minimizer : %d ",j-startIdx,data.minimizers[j]);
        //     }
        // }

        data.supermer_offs[insert_index] = 0;
        insert_index ++;

        unsigned int prev_minimizer = data.minimizers[startIdx];
        unsigned int cur_minimizer = data.minimizers[startIdx];

        int consecutive = 1;

        for(int i = startIdx + 1;i < startIdx + n;i++){

            // update two varaibles
            prev_minimizer = cur_minimizer;
            cur_minimizer  = data.minimizers[i];
            
            if (cur_minimizer == 88888888) { // reaches the endxw
                data.supermer_offs[insert_index] = n - K_kmer + 1;
                insert_index ++;
                break;
            }
            
            if (cur_minimizer == prev_minimizer){  
                consecutive++;
            }else{
                data.supermer_offs[insert_index] = data.supermer_offs[insert_index-1] + consecutive;
                insert_index ++;
                consecutive = 1; // reset 
            }
            
        }

        // if (index == 1){
        //     for(int i = data.reads_offs[index]; i < data.reads_offs[index] + n; i++){
        //         printf("%d:%d ",i-data.reads_offs[index],data.supermer_offs[i]);
        //     }
        // }
        
        // return;

        index += nthread;
    }
}


void GenerateSupermer_GPU(vector<string> &reads, int K, int P, vector<string> &all_supermers, int NUM_BLOCKS_PER_GRID, int NUM_THREADS_PER_BLOCK) {


    T_GPU_data d_batch_data, h_batch_data;
    T_read_count cur_batch_size;

    CSR<char, int> csr_reads(true);
        
    // convert raw read from vector to csr:
    int name_idx_in_batch = 0;
    for(string read: reads) {
        csr_reads.append(read.c_str(), read.length(),name_idx_in_batch);
        name_idx_in_batch++;
    }
    cur_batch_size = reads.size();
    // malloc and memcpy from host to device
    h_batch_data.reads = new char[csr_reads.size()];
    h_batch_data.reads_offs = csr_reads.get_raw_offs();
    h_batch_data.read_len = new T_read_len[csr_reads.items()];
    h_batch_data.minimizers = new T_minimizer[csr_reads.size()];
    h_batch_data.supermer_offs = new T_read_len[csr_reads.size()];
    for (int i=0; i<name_idx_in_batch-1; i++) {
        h_batch_data.read_len[i] = csr_reads.get_raw_offs()[i+1] - csr_reads.get_raw_offs()[i];
    }
    h_batch_data.read_len[name_idx_in_batch-1] = csr_reads.size() - csr_reads.get_raw_offs()[name_idx_in_batch-1];

    // [Data H2D]
    // gpu malloc
    GPUErrChk(cudaMalloc((void**) &(d_batch_data.reads), sizeof(char) * csr_reads.size()));
    GPUErrChk(cudaMalloc((void**) &(d_batch_data.reads_offs), sizeof(size_t) * (csr_reads.items()+1)));
    d_batch_data.cur_batch_size = cur_batch_size;
    #ifdef DEBUG
    // cerr << "(GPU"<<GPU_ID<<"): batch_size = " << cur_batch_size << ", bases = " << csr_reads.size() << endl;
    // cerr << " <G" << GPU_ID << "> " << cur_batch_size << "|" << csr_reads.size();
    #endif
    GPUErrChk(cudaMalloc((void**) &(d_batch_data.read_len), sizeof(T_read_len) * csr_reads.items()));
    GPUErrChk(cudaMalloc((void**) &(d_batch_data.minimizers), sizeof(T_minimizer) * csr_reads.size()));
    GPUErrChk(cudaMalloc((void**) &(d_batch_data.supermer_offs), sizeof(T_read_len) * csr_reads.size()));
    // memcpy Host -> Device
    GPUErrChk(cudaMemcpy(d_batch_data.reads, csr_reads.get_raw_data(), sizeof(char) * csr_reads.size(), cudaMemcpyHostToDevice));
    GPUErrChk(cudaMemcpy(d_batch_data.read_len, h_batch_data.read_len, sizeof(T_read_len) * (csr_reads.items()), cudaMemcpyHostToDevice));
    GPUErrChk(cudaMemcpy(d_batch_data.reads_offs, csr_reads.get_raw_offs(), sizeof(size_t) * (csr_reads.items()+1), cudaMemcpyHostToDevice));

    // [Computing]
	// update1:
	float kernel_time;
	cudaEvent_t cuda_start, cuda_end;
	cudaEventCreate(&cuda_start);
	cudaEventCreate(&cuda_end);
	cudaEventRecord(cuda_start);
	
    GPU_GenMinimizer<<<NUM_BLOCKS_PER_GRID, NUM_THREADS_PER_BLOCK/*, 0, cuda_stream*/>>>(d_batch_data, K, P);
    GPU_GenSKM<<<NUM_BLOCKS_PER_GRID, NUM_THREADS_PER_BLOCK/*, 0, cuda_stream*/>>>(d_batch_data, K, P);

	// update1:
    cudaDeviceSynchronize();
	cudaEventRecord(cuda_end);
	cudaEventSynchronize(cuda_start);
	cudaEventSynchronize(cuda_end);
	cudaEventElapsedTime(&kernel_time, cuda_start, cuda_end);
	fprintf(stderr, "Driver Time: %.9lf s\n", kernel_time / pow(10, 3));

    // [Data D2H]
    
    GPUErrChk(cudaMemcpy(h_batch_data.reads, d_batch_data.reads, sizeof(char) * (csr_reads.size()), cudaMemcpyDeviceToHost));
    GPUErrChk(cudaMemcpy(h_batch_data.read_len, d_batch_data.read_len, sizeof(T_read_len) * (csr_reads.items()), cudaMemcpyDeviceToHost));
    GPUErrChk(cudaMemcpy(h_batch_data.minimizers, d_batch_data.minimizers, sizeof(T_minimizer) * (csr_reads.size()), cudaMemcpyDeviceToHost));
    GPUErrChk(cudaMemcpy(h_batch_data.supermer_offs, d_batch_data.supermer_offs, sizeof(T_read_len) * (csr_reads.size()), cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();
    // GPUErrChk(cudaStreamSynchronize(cuda_stream));
    GPUErrChk(cudaFree(d_batch_data.reads));
    GPUErrChk(cudaFree(d_batch_data.reads_offs));
    GPUErrChk(cudaFree(d_batch_data.minimizers));
    GPUErrChk(cudaFree(d_batch_data.read_len));
    GPUErrChk(cudaFree(d_batch_data.supermer_offs));
    cudaDeviceSynchronize();
    
    for (int i=0; i<name_idx_in_batch; i++) {
        int skm_idx = 1;
        T_read_len *skm_offs = &h_batch_data.supermer_offs[h_batch_data.reads_offs[i]];
        vector<string>  supermers_local_all;
        T_read_len len = h_batch_data.read_len[i];
        while (*(skm_offs+skm_idx-1) != len- K +1) {
            int skm_len = *(skm_offs+skm_idx) - *(skm_offs+skm_idx-1) + K-1;
            char* t = (char *) malloc(skm_len+1);
            memcpy(t, h_batch_data.reads +h_batch_data.reads_offs[i]+ *(skm_offs+skm_idx-1), skm_len);
            t[skm_len] = '\0';
            supermers_local_all.push_back(t);
            skm_idx++;
        }

        all_supermers.insert(all_supermers.end(), supermers_local_all.begin(), supermers_local_all.end());
    }
    

    delete [] h_batch_data.reads;
    delete [] h_batch_data.read_len;
    delete [] h_batch_data.minimizers;
    delete [] h_batch_data.supermer_offs;


    return;
}