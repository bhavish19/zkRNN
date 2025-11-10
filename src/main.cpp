#include <stdio.h>
#include <vector>
#include <polynomial.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <sys/resource.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <string>
#include <time.h>
#include "mimc.h"
#include "quantization.h"
#include "GKR.h"
#include <chrono>
#include "utils.hpp"
#include "pol_verifier.h"
#include "config_pc.hpp"
#include "poly_commit.h"
#include "logup.hpp"
#include "bench.hpp"
#include "verifier.h"
#include "proof_serialization.h"
#include "witness_snapshot.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include<unistd.h>
using namespace std;

extern void   reset_gkr_verifier_time();
extern double get_gkr_verifier_time();

extern int partitions;

double g_verify_proof_time_last = 0.0;
double g_verify_proof_time_acc  = 0.0;

double g_commit_general = 0.0;
double g_commit_logup   = 0.0;
double g_commit_aggr    = 0.0;

int PC_scheme,Commitment_hash;
int levels;

double logup_acc_commit_time = 0.0;
double logup_acc_prove_time = 0.0;
double logup_acc_total_time = 0.0;

double forward_gkr_time = 0.0;
double aggregation_prove_time = 0.0;
double aggregation_time = 0.0;

double total_prove = 0.0;
double total_foward = 0.0;
double total_backward = 0.0;
double total_update = 0.0;
double total_logup = 0.0;
double total_logup_prove = 0.0;
double total_pogd = 0.0;
double total_aggregation_time = 0.0;
double total_aggregation_prove_time = 0.0;

unsigned long int mul_counter = 0;
double Forward_propagation_rnn;
double Backward_propagation_rnn;
double parameter_update_rnn;


vector<int> predicates_size;
vector<struct proof> Transcript;
vector<struct proof> Opening_proof;
vector<F> SHA;
vector<F> H;
extern int in_size;
int arity = 3;
vector<F> x_transcript,y_transcript;
F current_randomness;

static bool load_inputs_from_file(vector<vector<F>>& X,
                                  const std::string& path,
                                  int T,
                                  int n) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "[infer] Unable to open input file: %s\n", path.c_str());
        return false;
    }

    const size_t expected = static_cast<size_t>(T) * static_cast<size_t>(n);
    std::vector<float> flat;
    flat.reserve(expected);

    double value;
    while (flat.size() < expected && (in >> value)) {
        flat.push_back(static_cast<float>(value));
    }

    if (flat.size() != expected) {
        std::fprintf(stderr,
                     "[infer] Expected %zu input scalars but read %zu from %s\n",
                     expected,
                     flat.size(),
                     path.c_str());
        return false;
    }

    X.assign(T, std::vector<F>(n, F_ZERO));
    size_t idx = 0;
    for (int t = 0; t < T; ++t) {
        for (int j = 0; j < n; ++j) {
            X[t][j] = quantize(flat[idx++]);
        }
    }

    if (in >> value) {
        std::fprintf(stderr,
                     "[infer] Warning: extra values detected in %s; ignoring the tail\n",
                     path.c_str());
    }

    return true;
}

static bool load_model_from_file(rnn_layer& net,
                                 const std::string& path,
                                 int n,
                                 int m,
                                 int k) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "[infer] Unable to open model file: %s\n", path.c_str());
        return false;
    }

    const size_t expected = static_cast<size_t>(n) * static_cast<size_t>(m)
                          + static_cast<size_t>(m) * static_cast<size_t>(m)
                          + static_cast<size_t>(k) * static_cast<size_t>(m)
                          + static_cast<size_t>(m)
                          + static_cast<size_t>(k);

    std::vector<float> values;
    values.reserve(expected);

    double value;
    while (values.size() < expected && (in >> value)) {
        values.push_back(static_cast<float>(value));
    }

    if (values.size() != expected) {
        std::fprintf(stderr,
                     "[infer] Expected %zu model scalars but read %zu from %s\n",
                     expected,
                     values.size(),
                     path.c_str());
        return false;
    }

    size_t idx = 0;

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            net.W_x[i][j] = quantize(values[idx++]);
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            net.W_h[i][j] = quantize(values[idx++]);
        }
    }

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < m; ++j) {
            net.W_y[i][j] = quantize(values[idx++]);
        }
    }

    for (int i = 0; i < m; ++i) {
        net.b1[i] = quantize(values[idx++]);
    }

    for (int i = 0; i < k; ++i) {
        net.b2[i] = quantize(values[idx++]);
    }

    if (in >> value) {
        std::fprintf(stderr,
                     "[infer] Warning: extra values detected in %s; ignoring the tail\n",
                     path.c_str());
    }

    return true;
}

double PeakMemGB() {
  unsigned long long bytes = 0ULL;

  // 1) 优先从 /proc/self/status 读取 VmHWM（历史峰值RSS，单位 kB）
  if (std::ifstream in{"/proc/self/status"}) {
    std::string line;
    while (std::getline(in, line)) {
      if (line.rfind("VmHWM:", 0) == 0) {
        const char* s = line.c_str();
        while (*s && (*s < '0' || *s > '9')) ++s;        // 跳过到数字
        unsigned long long kb = std::strtoull(s, nullptr, 10);
        bytes = kb * 1024ULL;
        break;
      }
    }
  }

  // 2) 兜底：getrusage(RUSAGE_SELF).ru_maxrss（Linux 下单位 kB）
  if (bytes == 0ULL) {
    rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
      bytes = static_cast<unsigned long long>(ru.ru_maxrss) * 1024ULL;
    }
  }

  double gb = static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
  return gb; // 保留两位小数
}


void init_SHA(){
	for(int i = 0; i < 64; i++){
		SHA.push_back(F(random()%(1ULL<<32)));
	}
	for(int i = 0; i < 8; i++){
		H.push_back(F(random()%(1ULL<<32)));
	}
}


vector<vector<F>> prepare_matrixes(vector<vector<F>> &M1, vector<vector<F>> &M2, vector<F> r1, vector<F> r2){
	vector<vector<F>> V;
	
	V.push_back(prepare_matrix(transpose(M1),r1));

	V.push_back(prepare_matrix(transpose(M2),r2));
	
	return V;

}

struct proof generate_4product_sumcheck_proof(vector<F> &v1, vector<F> &v2,F previous_r){
	struct proof Pr;
	Pr.type = HASH_SUMCHECK;
	//vector<F> r = generate_randomness(int(log2(v1.size())));
	int rounds = int(log2(v1.size()));
	vector<quadruple_poly> p;
	F rand = previous_r;
	vector<F> r;
	for(int i = 0; i < rounds; i++){
		quadruple_poly poly = quadruple_poly(F_ZERO,F_ZERO,F_ZERO,F_ZERO,F_ZERO);
		cubic_poly temp_poly = cubic_poly(F_ZERO,F_ZERO,F_ZERO,F_ZERO);
		
		linear_poly l1,l4;
		
		int L = 1 << (rounds -1- i);
		for(int j = 0; j < L; j++){
			l1 = linear_poly(v1[2*j+1] - v1[2*j],v1[2*j]);
			//q1 = quadratic_poly()
			l4 = linear_poly(v2[2*j+1] - v2[2*j],v2[2*j]);
			//l3 = linear_poly(v3[2*j+1] - v3[2*j],v3[2*j]);
			temp_poly = l1*l1*l1;
			poly = poly + temp_poly*l4;

			v1[j] = v1[2*j] + rand*(v1[2*j+1]-v1[2*j]);
			v2[j] = v2[2*j] + rand*(v2[2*j+1]-v2[2*j]);
		}
		r.push_back(rand);
		vector<F> input;
		mimc_hash(poly.a,current_randomness);
		mimc_hash(poly.b,current_randomness);
		mimc_hash(poly.c,current_randomness);
		mimc_hash(poly.d,current_randomness);
		mimc_hash(poly.e,current_randomness);
		rand = current_randomness;
		p.push_back(poly);
	}
	Pr.quad_poly = p;
	Pr.randomness.push_back(r);
	mimc_hash(v1[0],current_randomness);
	mimc_hash(v2[0],current_randomness);
	
	Pr.vr.push_back(v1[0]);
	Pr.vr.push_back(v2[0]);
	
	return Pr;
}


void extend_input(vector<F> input, vector<F> extended_input, int partitions){
	vector<F> buff;
	for(int i = 0; i < input.size()/2; i++){
		buff = mimc_hash_segments(input[2*i],input[2*i+1]);
		for(int j = 0; j < partitions; j++){
			extended_input[j*input.size()/2 + i] = buff[j];
		}
	}
}


vector<proof> mimc_sumcheck(vector<F> input){
	proof P;
	int rounds = 80/partitions;
	pad_vector(input);
	vector<F> extended_input(input.size()*partitions/2);
	extend_input(input,extended_input,partitions);
	vector<vector<F>> data(rounds);
	vector<vector<F>> C_M(partitions);
	vector<F> C = get_parameters();
	
	clock_t start = clock();
	int counter = 0;
	for(int i = 0; i < partitions; i++){
		C_M[i].resize(rounds);
		for(int j = 0; j < rounds; j++){
			C_M[i][j] = C[counter];
			counter++;
		}
	}
	
	data[0] = extended_input;
	for(int i = 1; i < rounds; i++){
		data[i].resize(extended_input.size());
		for(int j = 0; j < input.size()/2; j++){
			for(int k = 0; k < partitions; k++){
				data[i][ k*input.size()/2 + j] = (data[i-1][k*input.size()/2 + j] + C_M[k][i])*(data[i-1][k*input.size()/2 + j] + C_M[k][i])*(data[i-1][k*input.size()/2 + j] + C_M[k][i]);  
			}
		}
	}
	
	
	vector<F> x = generate_randomness((int)log2(extended_input.size()),current_randomness);
	F sum = evaluate_vector(data[rounds-1],x);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	vector<F> beta;
	vector<F> v(data[0].size());
	vector<proof> proofs;

	for(int i = rounds-2; i >= 0; i--){
		start = clock();
		precompute_beta(x,beta);

		vector<F> temp_c;
		for(int j = 0; j < partitions; j++){
			temp_c.push_back(C_M[j][i+1]);
		}
		for(int j = 0; j < input.size()/2; j++){
			for( int k = 0; k < partitions; k++){
				v[k*input.size()/2 + j] = data[i][k*input.size()/2 + j] + C_M[k][i+1];
			}
		}
		total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;
		
		start = clock();
		P = generate_4product_sumcheck_proof(v, beta,current_randomness);
		total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

		proofs.push_back(P);
		if(P.quad_poly[0].eval(F(0))  +P.quad_poly[0].eval(F(1)) != sum){
			printf("Error %d\n",i);
		}
		// check correctness
		x = P.randomness[0];
		
		vector<F> x_c;
		for(int j = x.size()-(int)log2(partitions); j < x.size(); j++ ){
			x_c.push_back(x[j]);
		}

		start = clock();
		sum = P.vr[0] - evaluate_vector(temp_c,x_c);
		total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;
	}
	return proofs;

}


struct proof generate_3product_sumcheck_proof(vector<F> &v1, vector<F> &v2, vector<F> &v3,F previous_r){
	struct proof Pr;
	vector<F> r_sc;
	//vector<F> r = generate_randomness(int(log2(v1.size())));
	int rounds = int(log2(v1.size()));
	vector<cubic_poly> p;
	F rand = previous_r;

	for(int i = 0; i < rounds; i++){
		cubic_poly poly = cubic_poly(F_ZERO,F_ZERO,F_ZERO,F_ZERO);
		cubic_poly temp_poly = cubic_poly(F_ZERO,F_ZERO,F_ZERO,F_ZERO);
		linear_poly l1,l2,l3;
		
		int L = 1 << (rounds -1- i);
		for(int j = 0; j < L; j++){
			if(!(v2[2*j] == F(0) && v2[2*j+1] == F(0))){
				l1 = linear_poly(v1[2*j+1] - v1[2*j],v1[2*j]);
				//q1 = quadratic_poly()
				l2 = linear_poly(v2[2*j+1] - v2[2*j],v2[2*j]);
				l3 = linear_poly(v3[2*j+1] - v3[2*j],v3[2*j]);
				temp_poly = l1*l2*l3;
				poly = poly + temp_poly;

			}
			v1[j] = v1[2*j] + rand*(v1[2*j+1]-v1[2*j]);
			if(v2[2*j] == F(0) && v2[2*j+1] == F(0)){
				v2[j] = F(0);
				v3[j] = F(1);
			}else{
				v2[j] = v2[2*j] + rand*(v2[2*j+1]-v2[2*j]);
				v3[j] = F(1)-v2[j];//(1-r[i])*v3[2*j] + r[i]*v3[2*j+1];	

			}
		}
		r_sc.push_back(rand);
		vector<F> input;
		mimc_hash(poly.a,current_randomness);
		mimc_hash(poly.b,current_randomness);
		mimc_hash(poly.c,current_randomness);
		mimc_hash(poly.d,current_randomness);
		rand = current_randomness;
		p.push_back(poly);
	}
	Pr.c_poly = p;
	mimc_hash(v2[0],current_randomness);
	mimc_hash(v1[0],current_randomness);
	
	Pr.vr.push_back(v2[0]);
	Pr.vr.push_back(v1[0]);
	Pr.randomness.push_back(r_sc);
	Pr.sc_challenges.push_back(r_sc);


	return Pr;
}

struct proof generate_2product_sumcheck_proof(vector<F> v1, vector<F> v2, F previous_r){
	struct proof Pr;
	vector<F> r_sc;// = generate_randomness(int(log2(v1.size())));
	F rand = mimc_hash(previous_r,F(0));
	vector<quadratic_poly> p;
	int rounds = int(log2(v1.size()));
	for(int i = 0; i < rounds; i++){
		quadratic_poly poly = quadratic_poly(F_ZERO,F_ZERO,F_ZERO);
		quadratic_poly temp_poly = quadratic_poly(F_ZERO,F_ZERO,F_ZERO);
		linear_poly l1,l2;
		int L = 1 << (rounds - 1-i);
		for(int j = 0; j < L; j++){
			l1 = linear_poly(v1[2*j+1] - v1[2*j],v1[2*j]);
			l2 = linear_poly(v2[2*j+1] - v2[2*j],v2[2*j]);
			temp_poly = l1*l2;
			poly = poly + temp_poly;
			v1[j] = v1[2*j] + rand*(v1[2*j+1]-v1[2*j]);
			v2[j] = v2[2*j] + rand*(v2[2*j+1]-v2[2*j]);
		
		}
		r_sc.push_back(rand);
		vector<F> input;
		mimc_hash(poly.a,current_randomness);
		mimc_hash(poly.b,current_randomness);
		mimc_hash(poly.c,current_randomness);
		
		//vector<vector<F>> temp = mimc_multihash3(input);
		//Pr.w_hashes.push_back(temp);
		//rand = temp[temp.size()-1][temp[0].size()-1];
		//rand = mimc_multihash(input);
		rand = current_randomness;
		p.push_back(poly);
	}
	mimc_hash(v1[0],current_randomness);
	mimc_hash(v2[0],current_randomness);
	
	Pr.vr.push_back(v1[0]);
	Pr.vr.push_back(v2[0]);
	Pr.q_poly = p;
	Pr.randomness.push_back(r_sc);
	Pr.sc_challenges.push_back(r_sc);
	return Pr;
 }

struct proof _prove_matrix2matrix(vector<vector<F>> M1, vector<vector<F>> M2, vector<F> r_eval, F previous_sum){
	struct proof Pr;
	int r1_len = (int)log2(M1.size());
	int r2_len = (int)log2(M2.size());
	assert((int)r_eval.size() == r1_len + r2_len);

	vector<F> r2(r_eval.begin(), r_eval.begin()+ r2_len); 
	vector<F> r1(r_eval.begin()+r2_len, r_eval.end());

	
	clock_t start = clock();
	vector<vector<F>> V = prepare_matrixes(M1,M2,r1,r2);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	start = clock();
	if(V[0].size() != 1){
		Pr = generate_2product_sumcheck_proof(V[0],V[1],r_eval.back());
		Pr.randomness.push_back(r1);
		Pr.randomness.push_back(r2);

		if(previous_sum*F(1ULL<<Q) != (Pr.q_poly[0].eval(0) +Pr.q_poly[0].eval(1)) ){
			printf("Error in Matrix2Matrix multiplication\n");
			exit(-1);
		}
		Pr.type = MATMUL_PROOF;
	}else{
		
		Pr.randomness.push_back(r1);
		Pr.randomness.push_back(r2);
		if(previous_sum*F(1ULL<<Q) != V[0][0]*V[1][0]){
			printf("Error in Matrix2Matrix multiplication\n");
			exit(-1);
		}
		Pr.type = -1;
	}
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;
	return Pr;
}

vector<vector<F>> generate_bit_matrix(vector<F> bits,int domain){
	vector<vector<F>> M;
	int elements = bits.size()/domain; 
	M.resize(domain);
	for(int i = 0; i < M.size(); i++){
		for(int j = 0; j < elements; j++){
			M[i].push_back(bits[j*domain+i]);
		}
	}
	return M;
}

struct proof _prove_bit_decomposition(vector<F> bits, vector<F> r1, F previous_sum, int domain){
	//int domain = 256;
	vector<F> powers;
	powers.push_back(F(1));
	for(int i = 1; i < domain; i++){
		powers.push_back(F(2)*powers[i-1]);
	}

	//check_integrity(bits,num,powers);


	clock_t start,end;
			
	start = clock(); 
	vector<vector<F>> M = generate_bit_matrix(bits,domain);
	vector<F> v1 = prepare_matrix(M,r1);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	start = clock(); 
	struct proof Pr1 = generate_2product_sumcheck_proof(v1,powers,r1[r1.size()-1]);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	if(previous_sum != Pr1.q_poly[0].eval(0) + Pr1.q_poly[0].eval(1)){
		printf("Error in bit_decomposition\n");
		exit(-1);
	}
	vector<F> r2 = generate_randomness(int(log2(bits.size())),r1[r1.size()-1]);
	vector<F> beta;
	
	start = clock(); 
	precompute_beta(r2,beta);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	vector<F> inv_bits;
	for(int i  = 0; i < bits.size(); i++){
		inv_bits.push_back(F(1) - bits[i]);
	}

	start = clock(); 
	struct proof Pr2 = generate_3product_sumcheck_proof(beta,bits,inv_bits,r2[r2.size()-1]);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	struct proof Pr;
	Pr.randomness.push_back(r1);
	Pr.randomness.push_back(r2);
	Pr.randomness.push_back(Pr1.randomness[0]);
	Pr.randomness.push_back(Pr2.randomness[0]);
	Pr.vr.insert(Pr.vr.end(),Pr1.vr.begin(),Pr1.vr.end());
	Pr.vr.insert(Pr.vr.end(),Pr2.vr.begin(),Pr2.vr.end());
	Pr.q_poly = Pr1.q_poly;
	Pr.c_poly = Pr2.c_poly;
	Pr.type = RANGE_PROOF;
	Pr.w_hashes.insert(Pr.w_hashes.end(),Pr1.w_hashes.begin(),Pr1.w_hashes.end());
	Pr.w_hashes.insert(Pr.w_hashes.end(),Pr2.w_hashes.begin(),Pr2.w_hashes.end());
	
	return Pr;

}

struct proof prove_add(vector<F> X, vector<F> Y, const vector<F> &r_eval, F &previous_sum){
	int n = X.size();
	assert((int)Y.size() == n);
	clock_t start,end;
	start = clock();
	F sumX = evaluate_vector(X,r_eval);
	F sumY = evaluate_vector(Y,r_eval);
	end = clock();
	total_prove += (float)(end-start)/(float)CLOCKS_PER_SEC;
	if (previous_sum != sumX + sumY) {
        printf("Error in prove_add\n");
        exit(-1);
    }

	struct proof P;
	P.type = ADD_PROOF;
	P.in1 = sumX;
	P.in2 = sumY;
	P.out_eval = previous_sum;
	P.randomness.push_back(r_eval);

	return P;
}

struct proof prove_matrix_add(const vector<vector<F>>& A,
                              const vector<vector<F>>& B,
                              const vector<F>& r_eval,
                              const F& previous_sum) {
    int T = (int)A.size();
    int K = (int)A[0].size();
    assert((int)B.size() == T && (int)B[0].size() == K);

    vector<F> flatA = convert2vector(A);
    vector<F> flatB = convert2vector(B);

    clock_t start = clock();
    F sumA = evaluate_vector(flatA, r_eval);
    F sumB = evaluate_vector(flatB, r_eval);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    if (sumA + sumB != previous_sum) {
        printf("Error in prove_matrix_add\n");
        exit(-1);
    }

    struct proof P;
    P.type = ADD_PROOF;
    P.in1 = sumA;
    P.in2 = sumB;
    P.out_eval = previous_sum;
    P.randomness.push_back(r_eval);
    return P;
}

struct proof prove_update_vector_add(const vector<F>& b, const vector<F>& db, const vector<F>& bnext, const F& lr){
	struct proof P;
	vector<F> r;
	const F SCALE = F(1ULL<<Q);
	const F alpha = lr/SCALE;

	r = generate_randomness((int)log2(b.size()),F(0));

	clock_t start = clock();
	auto sum_b = evaluate_vector(b,r);
	auto sum_db = evaluate_vector(db,r);
	auto sum_b_next = evaluate_vector(bnext,r);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;
	if(sum_b_next != sum_b - alpha * sum_db){
		printf("Update linear relation failed.\n");
        exit(-1);
	}

	P.type = ADD_PROOF;
	P.in1 = sum_b;
	P.in2 = F_ZERO - alpha * sum_db;
	P.out_eval = sum_b_next;
	P.eval_point= r;

	return P;
}

struct proof prove_update_matrix_add(const vector<vector<F>>& W, const vector<vector<F>>& dW, const vector<vector<F>>& Wnext, const F& lr){
	struct proof P;
	vector<F> r;
	const F SCALE = F(1ULL<<Q);
	const F alpha = lr/SCALE;

	auto flat_W_next = convert2vector(Wnext);
	auto flat_dW = convert2vector(dW);
	auto flat_W = convert2vector(W);

	r = generate_randomness((int)log2(flat_W.size()), F(0));

	clock_t start = clock();
	auto sum_W_next = evaluate_vector(flat_W_next, r);
	auto sum_dW = evaluate_vector(flat_dW, r);
	auto sum_W = evaluate_vector(flat_W, r);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

	if(sum_W_next != sum_W - alpha * sum_dW){
		printf("Update linear relation failed.\n");
        exit(-1);
	}

	P.type = ADD_PROOF;
	P.in1 = sum_W;
	P.in2 = F_ZERO - alpha * sum_dW;
	P.out_eval = sum_W_next;
	P.eval_point= r;

	return P;
}

void prove_rnn_update(struct rnn_layer net){
	const double proof0   = total_prove;
	struct proof P;
	P = prove_update_matrix_add(net.W_h, net.bwd.dW_h, net.W_h_next, net.learning_rate);
	Transcript.push_back(P);
	P = prove_update_matrix_add(net.W_x, net.bwd.dW_x, net.W_x_next, net.learning_rate);
	Transcript.push_back(P);
	P = prove_update_matrix_add(net.W_y, net.bwd.dW_y, net.W_y_next, net.learning_rate);
	Transcript.push_back(P);
	P  = prove_update_vector_add(net.b1, net.bwd.db1, net.b1_next,net.learning_rate);
	Transcript.push_back(P);
	P  = prove_update_vector_add(net.b2, net.bwd.db2, net.b2_next,net.learning_rate);
	Transcript.push_back(P);
	const double  proof1         = total_prove     - proof0;
	

	parameter_update_rnn = proof1;

	total_update += parameter_update_rnn;
}

void prove_rnn_forward(struct rnn_layer net, vector<F> h_init){
	const double proof0   = total_prove;

	struct proof P;
	vector<F> r;
	F previous_sum;
	int T = net.seq_len;
	int m = net.hidden_size;
	int k = net.output_size;
	logup_reset_all_accumulators();
	vector<F> prev_h = h_init;
	//softmax lookup
	P = prove_logup(net.fwd.exp_softmax, 15, 32);
	Transcript.push_back(P);

	auto flat_z = convert2vector(net.fwd.z);
	r = generate_randomness((int)log2(k*T), F(0));

	clock_t start = clock();
	previous_sum = evaluate_vector(flat_z,r);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;

	vector<vector<F>> b2mat(T, vector<F>(k));
	for(int t = 0; t < T; ++t) {
    	b2mat[t] = net.b2;   
	}

	P = prove_matrix_add(net.fwd.Zy,b2mat,r,previous_sum);
	Transcript.push_back(P);

	previous_sum = P.in1;
	P = _prove_matrix2matrix(net.fwd.h,net.W_y,r,previous_sum);		
	Transcript.push_back(P);

	//tanh Lookup
	P = prove_logup(net.fwd.exp_tanh, 15,32);
	Transcript.push_back(P);

	auto flat_a = convert2vector(net.fwd.a);
	r = generate_randomness((int)log2(T*m),F(0));
	
	start = clock();
	previous_sum = evaluate_vector(flat_a,r);
	total_prove += (float)(clock()-start)/(float)CLOCKS_PER_SEC;


	vector<vector<F>> b1mat(T, vector<F>(m));
	for(int t = 0; t < T; ++t) {
    	b1mat[t] = net.b1;   
	}

	P = prove_matrix_add(net.fwd.WhWx,b1mat,r, previous_sum);
	Transcript.push_back(P);

	previous_sum = P.in1;

	P = prove_matrix_add(net.fwd.Wh, net.fwd.Wx, r, previous_sum);
	Transcript.push_back(P);

	F sum_Wh = P.in1;
	F sum_Wx = P.in2;

	P = _prove_matrix2matrix(net.fwd.x, net.W_x,r, sum_Wx);
    Transcript.push_back(P);
		
	vector<vector<F>> Hprev(T, vector<F>(m));
	for(int t = 0; t < T; ++t){
		if(t > 0) {
        	Hprev[t] = net.fwd.h[t-1];   // 用前一步 t-1 的 h
    	} else {
        	Hprev[t] = h_init;           // t=0 用初始化 h_init
    	}
	}
	P = _prove_matrix2matrix(Hprev,net.W_h,r,sum_Wh);
	Transcript.push_back(P);

	const double proof1 = total_prove - proof0;

	Forward_propagation_rnn = proof1;
	total_foward += Forward_propagation_rnn;
}

struct proof prove_rnn_db1(const rnn_layer& net){
	int T = net.seq_len;
	int m = net.hidden_size;

	auto M1 = transpose(net.bwd.da);
	vector<vector<F>> M2(1, vector<F>(T, F(1ULL<<Q)));

	const auto& flat = net.bwd.db1;
	int r1 = (int)log2(m);
	int R =  r1+1;
	auto r = generate_randomness(R, F(0));

	vector<F> r_flat(r.begin(), r.begin()+ r1);
	clock_t start = clock();
	F s = evaluate_vector(flat, r_flat);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

	return _prove_matrix2matrix(M1, M2, r, s);
}

struct proof prove_rnn_db2(const rnn_layer& net) {

    int T = net.seq_len;
    int k = net.output_size;

    auto M1 = transpose(net.bwd.delta_y);

    vector<vector<F>> M2(1, std::vector<F>(T, F(1ULL<<Q)));

    const auto& flat = net.bwd.db2;

    int r1 = (int)log2(k);
    int R  = r1 + 1;

    auto r = generate_randomness(R, F(0));
    vector<F> r_flat(r.begin(), r.begin()+ r1);
    
	clock_t start = clock();
	F s = evaluate_vector(flat, r_flat);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    return _prove_matrix2matrix(M1, M2, r, s);
}

struct proof prove_rnn_dWx(const rnn_layer& net) {
    int T = net.seq_len;
    int m = net.hidden_size;
    int n = net.input_size;

    auto flat_dw = convert2vector(net.bwd.dW_x);

    int bits = (int)log2(m*n);
    auto r = generate_randomness(bits, F(0));
	
	clock_t start = clock();
    F claimed_sum = evaluate_vector(flat_dw, r);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    auto M1 = transpose(net.bwd.da);   
    auto M2 = transpose(net.fwd.x); 

    return _prove_matrix2matrix(M1, M2, r, claimed_sum);
}

struct proof prove_rnn_dWh(const rnn_layer& net,const vector<F>& h_init) {
    int T = net.seq_len;
    int m = net.hidden_size;

    auto M1 = transpose(net.bwd.da);
    vector<vector<F>> prev_h_mat(T, vector<F>(m));
    for (int t = 0; t < T; ++t) {
        if (t > 0) {
            prev_h_mat[t] = net.fwd.h[t-1];
        } else {
            prev_h_mat[t] = h_init;
        }
    }
    auto M2 = transpose(prev_h_mat);

	vector<F> flat_dWh = convert2vector(net.bwd.dW_h);

    int bits = (int)log2(m*m);
    auto r = generate_randomness(bits, F(0));

	clock_t start = clock();
	F claimed_sum = evaluate_vector(flat_dWh, r);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    return _prove_matrix2matrix(M1, M2, r, claimed_sum);

}

struct proof prove_rnn_dWy(const rnn_layer& net) {

    int T = net.seq_len;
    int k = net.output_size;    
    int m = net.hidden_size;    
    auto flat_dw = convert2vector(net.bwd.dW_y);

    int bits = (int)log2(k*m);
    auto r_full = generate_randomness(bits, F(0));
	
	clock_t start = clock();
    F claimed_sum = evaluate_vector(flat_dw, r_full);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    auto M1 = transpose(net.bwd.delta_y);
    auto M2 = transpose(net.fwd.h);

    return _prove_matrix2matrix(M1, M2, r_full, claimed_sum);
}

void prove_rnn_dh(const rnn_layer& net, const vector<F>& r_eval){

	const int T = net.seq_len;
    const int m = net.hidden_size;
	auto flat_dh = convert2vector(net.bwd.dh);
	
	clock_t start = clock();
	F prev_sum = evaluate_vector(flat_dh,r_eval);
	total_prove += float(clock() - start) / CLOCKS_PER_SEC;

	vector<vector<F>> B(T, vector<F>(m, F_ZERO));
	for (int t = 0; t < T - 1; ++t) {
    	B[t] = net.bwd.WhT_da[t + 1];
	}
	proof Pr = prove_matrix_add(net.bwd.WyT_delta_y, B, r_eval, prev_sum);
	Transcript.push_back(Pr);
	F wy_sum = Pr.in1;
	F wh_sum = Pr.in2;

	Pr = _prove_matrix2matrix(net.bwd.delta_y, transpose(net.W_y),r_eval,wy_sum);
	Transcript.push_back(Pr);
	vector<vector<F>> DaShift(T, vector<F>(m, F_ZERO));
    for (int t = 0; t < T - 1; ++t) DaShift[t] = net.bwd.da[t + 1];

	Pr = _prove_matrix2matrix(DaShift, transpose(net.W_h),r_eval,wh_sum);
	Transcript.push_back(Pr);
}

void prove_tanh_backward(const rnn_layer& net, const vector<F>& r_eval){

    const F SCALE = F(1ULL << Q);
    const int T = net.seq_len;
    const int m = net.hidden_size;
    const int N = T * m;
    
    // flatten
    const auto flat_da = convert2vector(net.bwd.da);
    const auto flat_dh = convert2vector(net.bwd.dh);

    vector<F> flat_u;  
	flat_u.reserve(N);
	vector<F> flat_y2; 
	flat_y2.reserve(N);

    for (int t = 0; t < T; ++t) {
        const auto& tb = net.bwd.tanh_backprops[t];
        flat_u .insert(flat_u .end(), tb.one_minus_y2.begin(), tb.one_minus_y2.end());
		flat_y2.insert(flat_y2.end(), tb.y2.begin(), tb.y2.end());
    }

	clock_t start = clock();
    const F s_da  = evaluate_vector(flat_da, r_eval);
    const F target = s_da * SCALE;

    vector<F> beta;
	beta.reserve(N);
    precompute_beta(r_eval, beta);
    total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    vector<F> dh_beta(N);
    for (int i = 0; i < N; ++i) dh_beta[i] = flat_dh[i] * beta[i];
	
	start = clock();
    proof Pr = generate_2product_sumcheck_proof(dh_beta,flat_u,r_eval.back());
	total_prove += float(clock() - start) / CLOCKS_PER_SEC;

    Pr.type = MATMUL_PROOF;
    Pr.eval_point = r_eval;
    Transcript.push_back(Pr);

	start = clock();
    const F s_u = evaluate_vector(flat_u, r_eval);          

    vector<F> ones(N, SCALE);                      
    vector<F> neg_y2(N);
    for (int i = 0; i < N; ++i) neg_y2[i] = F_ZERO - flat_y2[i];
	total_prove += float(clock() - start) / CLOCKS_PER_SEC;
    
	Pr = prove_add(ones,neg_y2,r_eval,const_cast<F&>(s_u));
    Pr.type = ADD_PROOF;
	Pr.eval_point = r_eval;
    Transcript.push_back(Pr);
}

void prove_rnn_backward(struct rnn_layer net, vector<F> h_init){
	const double proof0   = total_prove;
	proof p_wh = prove_rnn_dWh(net, h_init);
    Transcript.push_back(p_wh);
    proof p_wx = prove_rnn_dWx(net);
    Transcript.push_back(p_wx);
    proof p_b1 = prove_rnn_db1(net);
    Transcript.push_back(p_b1);

	const int T = net.seq_len;
	const int m = net.hidden_size;

	vector<F> r_eval_tm = generate_randomness((int)log2(T*m), F_ZERO);

	prove_rnn_dh(net, r_eval_tm);
	prove_tanh_backward(net, r_eval_tm);

	proof p_b2 = prove_rnn_db2(net);
	Transcript.push_back(p_b2);
	proof p_wy = prove_rnn_dWy(net);
	Transcript.push_back(p_wy);

	//TODO delta_y

	
	const double proof1          = total_prove     - proof0;
	
	Backward_propagation_rnn = proof1;
	total_backward += Backward_propagation_rnn;
}

void check_input_rnn(int seq_len, int input_size){
	int T = seq_len;
	int n = input_size;
	clock_t start = clock();
	vector<vector<F>> X(T, vector<F>(n));
	for(int t = 0; t < T; t++){
		for (int j = 0; j < n; j++){
			X[t][j] = random();
		}
	}

	for(int t = 0; t < T; t++){
		for (int j = 0; j < n; j++){
			mimc_hash(X[t][j], current_randomness);
		}
		for(int j = 0; j <16; j++){
			mimc_hash(current_randomness, current_randomness);
		}
	}
	vector<F> input = x_transcript;
	x_transcript.clear();

	vector<proof> P = mimc_sumcheck(input);
	Transcript.insert(Transcript.end(),P.begin(),P.end());
	clock_t end = clock();
}

void get_witness_rnn(rnn_layer net, vector<F> &witness){
witness.clear();
    const auto& bt = net.fwd.exp_tanh;    // tanh 的 exp 查表批次
    const auto& bs = net.fwd.exp_softmax; // softmax 的 exp 查表批次

    // 预留空间（可选）
    witness.reserve(4 /*counts*/ + bt.Xq.size() + bt.Yq.size()
                                 + bs.Xq.size() + bs.Yq.size());

    // 写入计数，方便后续从 witness 中分段解析
    witness.push_back(F((unsigned long long)bt.Xq.size()));
    witness.push_back(F((unsigned long long)bt.Yq.size()));
    witness.push_back(F((unsigned long long)bs.Xq.size()));
    witness.push_back(F((unsigned long long)bs.Yq.size()));

    // 依次写入各批次的查询索引及对应取值（如果你保留了 Yq）
    witness.insert(witness.end(), bt.Xq.begin(), bt.Xq.end());
    witness.insert(witness.end(), bt.Yq.begin(), bt.Yq.end());
    witness.insert(witness.end(), bs.Xq.begin(), bs.Xq.end());
    witness.insert(witness.end(), bs.Yq.begin(), bs.Yq.end());
}


void get_model_rnn(rnn_layer net, vector<F> &model) {
    // W_x: [hidden_size][input_size]
    for(int i = 0; i < net.hidden_size; ++i) {
        model.insert(model.end(),
                     net.W_x[i].begin(),
                     net.W_x[i].end());
    }
    // W_h: [hidden_size][hidden_size]
    for(int i = 0; i < net.hidden_size; ++i) {
        model.insert(model.end(),
                     net.W_h[i].begin(),
                     net.W_h[i].end());
    }
    // W_y: [output_size][hidden_size]
    for(int i = 0; i < net.output_size; ++i) {
        model.insert(model.end(),
                     net.W_y[i].begin(),
                     net.W_y[i].end());
    }
    // b1: [hidden_size]
    model.insert(model.end(),
                 net.b1.begin(),
                 net.b1.end());
    // b2: [output_size]
    model.insert(model.end(),
                 net.b2.begin(),
                 net.b2.end());
}

namespace {

void CreateDirIfMissing(const std::string& dir) {
    if (dir.empty()) return;
#ifdef _WIN32
    _mkdir(dir.c_str());
#else
    mkdir(dir.c_str(), 0755);
#endif
}

void EnsureParentDir(const std::string& path) {
    auto sep = path.find_last_of("/\\");
    if (sep == std::string::npos) return;
    std::string dir = path.substr(0, sep);
    if (dir.empty()) return;

    std::string prefix;
    size_t pos = 0;

#ifdef _WIN32
    if (dir.size() >= 2 && dir[1] == ':') {
        prefix = dir.substr(0, 2);
        pos = 2;
    }
#endif

    for (; pos < dir.size(); ++pos) {
        char c = dir[pos];
        prefix.push_back(c);
        if (c == '/' || c == '\\') {
            CreateDirIfMissing(prefix);
        }
    }

    CreateDirIfMissing(dir);
}

std::string ResolveDumpPath() {
    if (const char* custom = std::getenv("BENCH_DUMP")) {
        if (custom[0] != '\0') return std::string(custom);
    }
    if (const char* bench_csv = std::getenv("BENCH_CSV")) {
        std::string path(bench_csv);
        const auto pos = path.find_last_of("/\\");
        if (pos != std::string::npos) {
            return path.substr(0, pos + 1) + "witness_dump.json";
        }
    }
    return std::string("../bench/witness_dump.json");
}

void DumpFieldVector(std::ostream& out, const vector<F>& vec) {
    out << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i) out << ", ";
        out << "\"" << FieldElementToString(vec[i]) << "\"";
    }
    out << "]";
}

void DumpFieldMatrix(std::ostream& out, const vector<vector<F>>& matrix, int indent) {
    out << "[\n";
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < indent; ++j) out.put(' ');
        DumpFieldVector(out, matrix[i]);
        if (i + 1 != matrix.size()) out << ",";
        out << "\n";
    }
    out << "]";
}

void DumpTrainingSnapshot(const rnn_layer& net,
                          const vector<vector<F>>& inputs,
                          const vector<F>& h0) {
    static bool dumped = false;
    if (dumped) return;

    std::string path = ResolveDumpPath();
    EnsureParentDir(path);

    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) {
        std::fprintf(stderr, "[dump] failed to open %s\n", path.c_str());
        return;
    }

    dumped = true;
    std::fprintf(stderr, "[dump] dumping snapshot to %s\n", path.c_str());
    out << "{\n";
    out << "  \"seq_len\": " << net.seq_len << ",\n";
    out << "  \"input_size\": " << net.input_size << ",\n";
    out << "  \"hidden_size\": " << net.hidden_size << ",\n";
    out << "  \"output_size\": " << net.output_size << ",\n";

    out << "  \"input_t0\": ";
    if (!inputs.empty()) {
        DumpFieldVector(out, inputs.front());
    } else {
        out << "[]";
    }
    out << ",\n";

    out << "  \"h_prev_t0\": ";
    DumpFieldVector(out, h0);
    out << ",\n";

    auto dump_optional_vector = [&](const char* label, const vector<vector<F>>& container) {
        out << "  \"" << label << "\": ";
        if (!container.empty()) {
            DumpFieldVector(out, container.front());
        } else {
            out << "[]";
        }
        out << ",\n";
    };

    dump_optional_vector("a_t0", net.fwd.a);
    dump_optional_vector("h_t0", net.fwd.h);
    dump_optional_vector("z_t0", net.fwd.z);
    dump_optional_vector("y_hat_t0", net.fwd.yHat);

    out << "  \"exp_tanh\": {\n";
    out << "    \"Xq\": ";
    DumpFieldVector(out, net.fwd.exp_tanh.Xq);
    out << ",\n";
    out << "    \"Yq\": ";
    DumpFieldVector(out, net.fwd.exp_tanh.Yq);
    out << "\n  },\n";

    out << "  \"exp_softmax\": {\n";
    out << "    \"Xq\": ";
    DumpFieldVector(out, net.fwd.exp_softmax.Xq);
    out << ",\n";
    out << "    \"Yq\": ";
    DumpFieldVector(out, net.fwd.exp_softmax.Yq);
    out << "\n  },\n";

    out << "  \"weights\": {\n";
    out << "    \"W_x\": ";
    DumpFieldMatrix(out, net.W_x, 6);
    out << ",\n";
    out << "    \"W_h\": ";
    DumpFieldMatrix(out, net.W_h, 6);
    out << ",\n";
    out << "    \"W_y\": ";
    DumpFieldMatrix(out, net.W_y, 6);
    out << ",\n";
    out << "    \"b1\": ";
    DumpFieldVector(out, net.b1);
    out << ",\n";
    out << "    \"b2\": ";
    DumpFieldVector(out, net.b2);
    out << "\n  }\n";

    out << "}\n";
    out.close();
    std::fprintf(stderr, "[dump] wrote snapshot to %s\n", path.c_str());
}

} // namespace

vector<F> rotate(vector<F> bits, int shift){
	vector<F> new_bits;
 	for(int i = shift; i < bits.size(); i++){
        new_bits.push_back(bits[i]);
    }
    for(int i = 0; i < shift; i++){
        new_bits.push_back(bits[i]);
    }
	return new_bits;
}

vector<F> shift(vector<F> bits, int shift){
	vector<F> new_bits;
	for(int i = shift; i < bits.size(); i++){
		new_bits.push_back(bits[i]);
	}
	for(int i = 0; i < bits.size()-shift; i++){
		new_bits.push_back(F(0));
	}
	return new_bits;
}

SHA_witness get_sha_witness(vector<F> words){
	SHA_witness gkr_input;
	vector<F> pow(32);pow[0] = F(1);pow[1] = F(2);
	for(int i = 2; i < 32; i++){
		pow[i] = pow[i-1]*pow[1];
	}
	vector<F> new_words(48);
	vector<vector<F>> words_bits(64);
	vector<F> buff;
	for(int i = 0; i < 16; i++){
		buff.push_back(words[i]);
		if(words[i].toint128() > 1ULL<<32){
			printf("Error in word %d\n",i );
			exit(-1);
		}
	
		words_bits[i] = prepare_bit_vector(buff, 32);
		buff.clear();
	}
	vector<F> quotients;

	for(int i = 0;i < 16; i++){
		quotients.push_back(F_ZERO);
	}
	for(int i = 16; i < 64; i++){
		F temp_word1 = F_ZERO,temp_word2 = F_ZERO;
		vector<F> w1 = rotate(words_bits[i-15],7);
		vector<F> w2 = rotate(words_bits[i-15],18);
		vector<F> w3 = shift(words_bits[i-15],3);
		for(int j = 0; j < 32; j++){
			F _t = w1[j]+w2[j] - F(2)*w1[j]*w2[j];
			temp_word1 = temp_word1 + pow[j]*(w3[j] + _t - F(2)*w3[j]*_t);
		}
		
		w1 = rotate(words_bits[i-2],17);
		w2 = rotate(words_bits[i-2],19);
		w3 = shift(words_bits[i-2],10);
		for(int j = 0; j < 32; j++){
			F _t = w1[j]+w2[j] - F(2)*w1[j]*w2[j];
			temp_word2 =temp_word2 + pow[j]*(w3[j] + _t - F(2)*w3[j]*_t);
		}
		F temp_word = temp_word1 + temp_word2 + words[i - 16] + words[i-7];
		quotients.push_back(F(temp_word.toint128()/((unsigned long long)1ULL<<32)));
		words.push_back(temp_word.toint128()%((unsigned long long)1ULL<<32));
		buff.push_back(words[i]);
		if(words[i].toint128() > 1ULL<<32){
			printf("Error in word %d\n",i);
			exit(-1);
		}
		words_bits[i] = prepare_bit_vector(buff, 32);
		if(words_bits[i].size() != 32){
			printf("Error in witness 0\n");
			exit(-1);
		}
		buff.clear();
	}

	
	vector<F> a(65),b(65),c(65),d(65),e(65),f(65),g(65),h(65);
    a[0] = H[0];b[0] = H[1];c[0] = H[2];d[0] = H[3];
	e[0] = H[4];f[0] = H[5];g[0] = H[6];h[0] = H[7];
	vector<vector<F>> a_bits(64),b_bits(64),c_bits(64),e_bits(64),f_bits(64),g_bits(64);
    vector<F> a_q,e_q;

    for(int i = 0; i < 64; i++){
		buff.clear();
		buff.push_back(a[i]);
		a_bits[i] = prepare_bit_vector(buff,32);
		buff.clear();
		F s0 = F(0);
		vector<F> w1 = rotate(a_bits[i],2);
		vector<F> w2 = rotate(a_bits[i],13);
		vector<F> w3 = rotate(a_bits[i],22);
		
		for(int j = 0; j < 32; j++){
			F _t = w1[j]+w2[j] - F(2)*w1[j]*w2[j];
			s0 = s0+ pow[j]*(w3[j] + _t - F(2)*w3[j]*_t);
		}
		buff.push_back(b[i]);
		b_bits[i] = prepare_bit_vector(buff,32);
		buff.clear();
		buff.push_back(c[i]);
		c_bits[i] = prepare_bit_vector(buff,32);
		buff.clear();	
		if(a_bits[i].size() != 32 || b_bits[i].size() != 32 || c_bits[i].size() != 32){
			printf("Error in witness 1\n");
			exit(-1);
		}
		F maj = F(0);
		for(int j = 0; j < 32; j++){
			F _t = a_bits[i][j]*b_bits[i][j] + b_bits[i][j]*c_bits[i][j] - F(2)*a_bits[i][j]*b_bits[i][j]*b_bits[i][j]*c_bits[i][j];
			maj = maj + pow[j]*(c_bits[i][j]*a_bits[i][j] + _t - F(2)*c_bits[i][j]*a_bits[i][j]*_t);
		}
		F t2 = maj + s0;
		buff.push_back(e[i]);
		e_bits[i] = prepare_bit_vector(buff,32);
		buff.clear();
		buff.push_back(f[i]);
		f_bits[i] = prepare_bit_vector(buff,32);
		buff.clear();
		buff.push_back(g[i]);
		g_bits[i] = prepare_bit_vector(buff,32);
		buff.clear();
		if(e_bits[i].size() != 32 || f_bits[i].size() != 32 || g_bits[i].size() != 32){
			printf("Error in witness 1\n");
			exit(-1);
		}
		w1 = rotate(e_bits[i],6);
		w2 = rotate(e_bits[i],11);
		w3 = rotate(e_bits[i],25);
		F s1 = F(0);
		for(int j = 0; j < 32; j++){
			F _t = w1[j]+w2[j] - F(2)*w1[j]*w2[j];
			s1 = s1 + pow[j]*(w3[j] + _t - F(2)*w3[j]*_t);
		}
		F ch = F(0);
		for(int j = 0; j < 32; j++){
			ch = ch+ pow[j]*(e_bits[i][j]*f_bits[i][j] + (F(1)-e_bits[i][j])*g_bits[i][i] - F(2)*e_bits[i][j]*f_bits[i][j]* (F(1)-e_bits[i][j])*g_bits[i][i]);
		}
		F t1 = ch + s1 + words[i] + h[i] + SHA[i];
		a_q.push_back((t1+t2).toint128()/(1ULL<<32));
		a[i+1] = (t1+t2).toint128()%((1ULL<<32));
		e_q.push_back((t1 + d[i]).toint128()/(1ULL<<32));
		e[i+1] = (t1 + d[i]).toint128()%(1ULL<<32);
		h[i+1] = g[i];
		g[i+1] = f[i];
		f[i+1] = e[i];
		d[i+1] = c[i];
		c[i+1] = b[i];
		b[i+1] = a[i];
	}
	


	gkr_input.aux.insert(gkr_input.aux.end(),words.begin(),words.end());
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),words_bits[i].begin(),words_bits[i].end());
	}
	gkr_input.numbers.insert(gkr_input.numbers.end(),words.begin(),words.end());
	gkr_input.aux.insert(gkr_input.aux.end(),quotients.begin(),quotients.end());
	//gkr_input.insert(gkr_input.end(),pow.begin(),pow.end());
	
	for(int i = 0; i < a.size(); i++){
		gkr_input.aux.push_back(a[i]);
		gkr_input.aux.push_back(b[i]);
		gkr_input.aux.push_back(c[i]);
		gkr_input.aux.push_back(d[i]);
		gkr_input.aux.push_back(e[i]);
		gkr_input.aux.push_back(f[i]);
		gkr_input.aux.push_back(g[i]);
		gkr_input.aux.push_back(h[i]);
	}
	for(int i = 0; i < a_q.size(); i++){
		gkr_input.aux.push_back(a_q[i]);
		gkr_input.aux.push_back(e_q[i]);
	}
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),a_bits[i].begin(),a_bits[i].end());
		gkr_input.numbers.push_back(a[i]);
	}
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),b_bits[i].begin(),b_bits[i].end());
		gkr_input.numbers.push_back(b[i]);
	}
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),c_bits[i].begin(),c_bits[i].end());
		gkr_input.numbers.push_back(c[i]);
	}
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),e_bits[i].begin(),e_bits[i].end());
		gkr_input.numbers.push_back(e[i]);
	}
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),f_bits[i].begin(),f_bits[i].end());
		gkr_input.numbers.push_back(f[i]);
	}
	for(int i = 0; i < 64; i++){
		gkr_input.bits.insert(gkr_input.bits.end(),g_bits[i].begin(),g_bits[i].end());
		gkr_input.numbers.push_back(g[i]);
	}

	//gkr_input.push_back(F(1));
	return gkr_input;
}


vector<F> prove_aggregation(aggregation_witness data, int level){
	vector<F> witness;
	vector<F> gkr_input,numbers;
	gkr_input.insert(gkr_input.end(),SHA.begin(),SHA.end());
	vector<F> pow(32);pow[0] = F(1);pow[1] = F(2);
	for(int i = 2; i < 32; i++){
		pow[i] = pow[i-1]*pow[1];
	}
	gkr_input.insert(gkr_input.end(),pow.begin(),pow.end());
	vector<F> bits;
	printf("Merkle tree size : %d, %d\n",data.merkle_proof.size(),data.merkle_proof.size()/(arity+1));
	for(int i = 0; i < data.merkle_proof.size()/(arity+1); i+=2){
		
		vector<F> words;
		
		for(int j = 0; j < 8; j++){
			words.push_back(F(data.merkle_proof[i][j]));
		}
		for(int j = 0; j < 8; j++){
			words.push_back(F(data.merkle_proof[i+1][j]));
		}
		SHA_witness temp = get_sha_witness(words);
		gkr_input.insert(gkr_input.end(),temp.aux.begin(),temp.aux.end()); 
		gkr_input.insert(gkr_input.end(),temp.bits.begin(),temp.bits.end()); 
		bits.insert(bits.end(),temp.bits.begin(),temp.bits.end());
		numbers.insert(numbers.end(),temp.numbers.begin(),temp.numbers.end());
	}
	witness = gkr_input;
	gkr_input.push_back(F(1));

	if(data.merkle_proof.size()/(2*(arity+1)) != 0){
		struct proof _P;
		
		_P = prove_sha256(gkr_input, generate_randomness(10,F(312)),data.merkle_proof.size()/(2*(arity+1)));
		for(int i = 0; i < arity+1; i++){
			Transcript.push_back(_P);
		}
		gkr_input.clear();
		for(int i = 0; i < data.merkle_proof.size()/(arity+1); i++){
			for(int j = 0; j < data.merkle_proof[i].size(); j++){
				gkr_input.push_back(F(data.merkle_proof[i][j]));
			}
		}
		for(int i = 0; i < data.output.size()/(arity+1); i++){
			for(int j = 0; j < data.output[i].size(); j++){
				gkr_input.push_back(data.output[i][j]);
			}
		}
		for(int i = 0; i < data.idx.size()/(arity+1); i++){
			gkr_input.push_back(data.idx[i]);	
		}
		
		for(int i = 0; i < data.root_sha.size()/(arity+1); i++){
			for(int j = 0; j < data.root_sha[i].size(); j++){
				gkr_input.push_back(data.root_sha[i][j]);
			}
		}
		witness.insert(witness.end(),gkr_input.begin(),gkr_input.end());
		
		gkr_input.push_back(data.a[0]);
		gkr_input.push_back(F(-1));
		
		
		gkr_input.push_back(F_ONE);
		gkr_input.push_back(F(1ULL<<32));

		_P = prove_merkle_proof_consistency(gkr_input,generate_randomness(10,F(9)),data.merkle_proof.size()/2,level,data.trees);

		for(int i = 0; i < arity+1; i++){
			Transcript.push_back(_P);
		}
		gkr_input.clear();

		pad_vector(bits);pad_vector(numbers);
		vector<F> r = generate_randomness((int)log2(numbers.size()),F(32));
		_P = _prove_bit_decomposition(bits,r,evaluate_vector(numbers,r),32);

		for(int i = 0; i  < arity+1; i++){
			Transcript.push_back(_P);
		}
	}
	vector<F> one_node_hashes;
	for(int i = 0; i < data.merkle_proof_f.size()/(arity+1); i++){
		one_node_hashes.push_back(data.merkle_proof_f[i]);
	}
	vector<proof> P = mimc_sumcheck(one_node_hashes);
	Transcript.insert(Transcript.end(),P.begin(),P.end());
	for(int i = 0; i < arity; i++){
		Transcript.insert(Transcript.end(),P.begin(),P.end());
	}
	return witness;
}

// Note that prove aggr needs to return a witness w.
vector<F> prove_aggr(vector<vector<vector<vector<F>>>> &matrixes,vector<vector<commitment>> &comm){
	// Get the randomness to aggregate polynomials
	vector<F> total_witness,temp_w;
	vector<proof> aggr_proof;
	vector<F> a = generate_randomness(2,current_randomness);
	aggregation_witness aggregation_data;
	printf("Start aggregation\n");
	for(int i = 0; i < matrixes.size(); i++){
		// first measure the aggregation time
		aggregation_time = 0;
		aggregation_data = (aggregate(matrixes[i], comm[i],  a,  levels));
		total_aggregation_time += aggregation_time;
		// measure the recursion time
		clock_t start = clock();
		aggregation_prove_time = 0;
		temp_w = prove_aggregation(aggregation_data,levels);
		aggregation_prove_time = float(clock() - start) / CLOCKS_PER_SEC;
		total_aggregation_prove_time += aggregation_prove_time;
		total_witness.insert(total_witness.end(),temp_w.begin(),temp_w.end());
	}
	return total_witness;
	
}

void adopt_next_and_drop_grads(rnn_layer& net) {
    // 1) 采纳 next -> current（零拷贝）
    net.W_x.swap(net.W_x_next);
    net.W_h.swap(net.W_h_next);
    net.W_y.swap(net.W_y_next);
    net.b1.swap(net.b1_next);
    net.b2.swap(net.b2_next);

    // 2) 释放 *_next 的容量
    clear_matrix(net.W_x_next);
    clear_matrix(net.W_h_next);
    clear_matrix(net.W_y_next);
    clear_arr(net.b1_next);
    clear_arr(net.b2_next);

    // 3) 清本轮梯度（下一轮会重算）
    clear_matrix(net.bwd.dW_x);  { std::vector<std::vector<F>>().swap(net.bwd.dW_x); }
    clear_matrix(net.bwd.dW_h);  { std::vector<std::vector<F>>().swap(net.bwd.dW_h); }
    clear_matrix(net.bwd.dW_y);  { std::vector<std::vector<F>>().swap(net.bwd.dW_y); }
    clear_arr(net.bwd.db1);      { std::vector<F>().swap(net.bwd.db1); }
    clear_arr(net.bwd.db2);      { std::vector<F>().swap(net.bwd.db2); }

    // 4) 清运行期缓存（前向/后向临时量、LogUp批等）
    clear_matrix(net.fwd.x);
    clear_matrix(net.fwd.a);
    clear_matrix(net.fwd.h);
    clear_matrix(net.fwd.z);
    clear_matrix(net.fwd.yHat);
    clear_matrix(net.fwd.Wh);
    clear_matrix(net.fwd.Wx);
    clear_matrix(net.fwd.WhWx);
    clear_matrix(net.fwd.Zy);
    vector<tanh_layer>().swap(net.fwd.tanhs);
    vector<softmax_layer>().swap(net.fwd.softmaxes);
    net.fwd.exp_tanh = ExpBatch{};
    net.fwd.exp_softmax = ExpBatch{};
	clear_matrix(net.bwd.delta_y);
    clear_matrix(net.bwd.dh);
    clear_matrix(net.bwd.da);
    clear_matrix(net.bwd.WyT_delta_y);
    clear_matrix(net.bwd.WhT_da);
    vector<tanh_layer_backprop>().swap(net.bwd.tanh_backprops);
}

vector<F> init_row(int dim)
{
    vector<F> row(dim);
    for (int i = 0; i < dim; ++i)
    {
        row[i] = quantize((float)rand() / (float)RAND_MAX);
    }
    return row;
}

vector<vector<F>> init_input(int T, int n)
{
    vector<vector<F>> in(T);
    for (int t = 0; t < T; ++t)
    {
        in[t] = init_row(n);
    }
    return in;
} 

int main(int argc, char *argv[]){
    if (argc < 7) {
        std::fprintf(stderr,
                     "Usage: %s T n m k pc_type arity [--mode=train|infer] [--model=path] [--inputs=path]\n",
                     argv[0]);
        return 1;
    }

    init_hash();
    init_SHA();
    PC_scheme = 1;
    Commitment_hash = 1;
    levels = 2;
    
    int T = atoi(argv[1]); //seq_len
    int n = atoi(argv[2]);//input_size
    int m = atoi(argv[3]);//hidden_size
    int k = atoi(argv[4]); //output_size
    PC_scheme = atoi(argv[5]);
    arity = atoi(argv[6]);

    std::string mode = "train";
    std::string model_path;
    std::string input_path;
    std::string snapshot_path;

    for (int i = 7; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg.rfind("--mode=", 0) == 0) {
            mode = arg.substr(7);
        } else if (arg.rfind("--model=", 0) == 0) {
            model_path = arg.substr(8);
        } else if (arg.rfind("--inputs=", 0) == 0) {
            input_path = arg.substr(9);
        } else if (arg.rfind("--snapshot=", 0) == 0) {
            snapshot_path = arg.substr(11);
        } else {
            std::fprintf(stderr, "Unknown argument: %s\n", arg.c_str());
            return 1;
        }
    }

    WitnessSnapshot snapshot;
    bool use_snapshot = false;
    if (!snapshot_path.empty()) {
        try {
            snapshot = LoadWitnessSnapshot(snapshot_path);
            use_snapshot = true;
        } catch (const std::exception &e) {
            std::fprintf(stderr, "[snapshot] Failed to load %s: %s\n",
                         snapshot_path.c_str(), e.what());
            return 1;
        }

        auto adopt_dim = [](int current, int candidate,
                            const char *label) -> int {
            if (candidate <= 0)
                return current;
            if (current != candidate) {
                std::fprintf(stderr,
                             "[snapshot] Overriding %s from %d to %d\n",
                             label, current, candidate);
            }
            return candidate;
        };

        T = adopt_dim(T, snapshot.seq_len, "seq_len");
        n = adopt_dim(n, snapshot.input_size, "input_size");
        m = adopt_dim(m, snapshot.hidden_size, "hidden_size");
        k = adopt_dim(k, snapshot.output_size, "output_size");
    }

    const bool is_inference = (mode == "infer" || mode == "inference");

    vector<vector<F>> X_rnn;
    FieldVector h0;
    if (use_snapshot) {
        if (T <= 0 || n <= 0 || m <= 0) {
            std::fprintf(stderr,
                         "[snapshot] Invalid dimensions after loading snapshot: "
                         "T=%d n=%d m=%d\n",
                         T, n, m);
            return 1;
        }
        X_rnn.assign(static_cast<size_t>(T), vector<F>(n, F_ZERO));
        if (!snapshot.input_t0.empty()) {
            if (static_cast<int>(snapshot.input_t0.size()) != n) {
                std::fprintf(stderr,
                             "[snapshot] input_t0 length mismatch: expected %d "
                             "entries, got %zu\n",
                             n, snapshot.input_t0.size());
                return 1;
            }
            if (!X_rnn.empty()) {
                X_rnn[0] = snapshot.input_t0;
                for (int t = 1; t < T; ++t) {
                    X_rnn[t] = snapshot.input_t0;
                }
            }
        } else {
            std::fprintf(stderr,
                         "[snapshot] input_t0 missing; using zeros for inputs\n");
        }

        h0.assign(m, F_ZERO);
        if (!snapshot.h_prev_t0.empty()) {
            if (static_cast<int>(snapshot.h_prev_t0.size()) != m) {
                std::fprintf(stderr,
                             "[snapshot] h_prev_t0 length mismatch: expected %d "
                             "entries, got %zu\n",
                             m, snapshot.h_prev_t0.size());
                return 1;
            }
            h0 = snapshot.h_prev_t0;
        }
    } else {
        if (!input_path.empty()) {
            if (!load_inputs_from_file(X_rnn, input_path, T, n)) {
                return 1;
            }
        } else {
            X_rnn = init_input(T, n);
        }
        h0.assign(m, F_ZERO);
    }

    struct rnn_layer rnn_net = rnn_init(T,n,m,k);

    if (!model_path.empty() && use_snapshot) {
        std::fprintf(stderr,
                     "[snapshot] Warning: both --model and --snapshot supplied; "
                     "ignoring model file in favour of snapshot\n");
    }

    if (!use_snapshot && !model_path.empty()) {
        if (!load_model_from_file(rnn_net, model_path, n, m, k)) {
            return 1;
        }
    }

    if (use_snapshot) {
        auto validate_matrix = [](const std::vector<FieldVector> &mat,
                                  int rows,
                                  int cols,
                                  const char *label) -> bool {
            if (mat.empty())
                return true;
            if (static_cast<int>(mat.size()) != rows) {
                std::fprintf(stderr,
                             "[snapshot] %s row mismatch: expected %d, got %zu\n",
                             label, rows, mat.size());
                return false;
            }
            for (size_t i = 0; i < mat.size(); ++i) {
                if (static_cast<int>(mat[i].size()) != cols) {
                    std::fprintf(stderr,
                                 "[snapshot] %s column mismatch at row %zu: "
                                 "expected %d, got %zu\n",
                                 label, i, cols, mat[i].size());
                    return false;
                }
            }
            return true;
        };

        if (!validate_matrix(snapshot.W_x, m, n, "W_x") ||
            !validate_matrix(snapshot.W_h, m, m, "W_h") ||
            !validate_matrix(snapshot.W_y, k, m, "W_y")) {
            return 1;
        }

        if (!snapshot.W_x.empty()) rnn_net.W_x = snapshot.W_x;
        if (!snapshot.W_h.empty()) rnn_net.W_h = snapshot.W_h;
        if (!snapshot.W_y.empty()) rnn_net.W_y = snapshot.W_y;

        if (!snapshot.b1.empty()) {
            if (static_cast<int>(snapshot.b1.size()) != m) {
                std::fprintf(stderr,
                             "[snapshot] b1 length mismatch: expected %d, got %zu\n",
                             m, snapshot.b1.size());
                return 1;
            }
            rnn_net.b1 = snapshot.b1;
        }
        if (!snapshot.b2.empty()) {
            if (static_cast<int>(snapshot.b2.size()) != k) {
                std::fprintf(stderr,
                             "[snapshot] b2 length mismatch: expected %d, got %zu\n",
                             k, snapshot.b2.size());
                return 1;
            }
            rnn_net.b2 = snapshot.b2;
        }
    }

    if (is_inference) {
        Transcript.clear();
        x_transcript.clear();

        total_prove = 0.0;
        Forward_propagation_rnn = 0.0;
        total_foward = 0.0;
        logup_acc_total_time = 0.0;
        logup_acc_prove_time = 0.0;
        total_logup = 0.0;
        total_logup_prove = 0.0;

        rnn_net = rnn_forward(rnn_net, X_rnn, h0);

        std::vector<F> witness_rnn;
        get_witness_rnn(rnn_net, witness_rnn);

        std::vector<F> model_rnn;
        get_model_rnn(rnn_net, model_rnn);
        witness_rnn.insert(witness_rnn.end(), model_rnn.begin(), model_rnn.end());
        pad_vector(witness_rnn);

        std::vector<std::vector<F>> witness_matrix_rnn;
        commitment commitment_rnn;
        double commit_before = g_commit_general;
        poly_commit(witness_rnn, witness_matrix_rnn, commitment_rnn, levels);
        double commit_time = g_commit_general - commit_before;

        prove_rnn_forward(rnn_net, h0);

        double prover_time = Forward_propagation_rnn;
        double logup_total = logup_acc_total_time;
        double logup_prove = logup_acc_prove_time;

        float forward_proof_size = proof_size(Transcript);
        printf("Inference proof size: %lf\n", forward_proof_size);
        printf("Inference prover time: %lf\n", prover_time);
        printf("Commitment time: %lf\n", commit_time);
        printf("LogUp time (total / prove): %lf / %lf\n", logup_total, logup_prove);

        reset_gkr_verifier_time();
        verify_proof(Transcript);
        double inference_verifier_time = get_gkr_verifier_time();
        printf("Inference verifier time: %lf\n", inference_verifier_time);

        double peak_mem = PeakMemGB();
        printf("Peak resident memory: %.3f GB\n", peak_mem);

        BenchmarkRow row;
        row.timestamp = now_iso8601();
        row.hostname = get_hostname();
        row.T = T; row.n = n; row.m = m; row.k = k;
        row.PC_scheme = PC_scheme;
        row.Commitment_hash = Commitment_hash;
        row.levels = levels;
        row.arity = 1;
        row.param_count = 1LL*n*m + 1LL*m*m + 1LL*k*m + m + k;
        row.data_prove_time = 0.0;
        row.avg_update = 0.0;
        row.avg_backward = 0.0;
        row.avg_forward = Forward_propagation_rnn;
        row.avg_logup = logup_acc_total_time;
        row.avg_logup_prove = logup_acc_prove_time;
        row.avg_pogd = 0.0;
        row.agg_recursion_total = 0.0;
        row.agg_recursion_amort = 0.0;
        row.aggr_time_total = 0.0;
        row.aggr_time_avg = 0.0;
        row.aggr_prove_total = 0.0;
        row.aggr_prove_avg = 0.0;
        row.commit_time_avg = commit_time;
        row.recursion_verifier_time = inference_verifier_time;
        row.verifier_time = inference_verifier_time;
        row.pogd_proof_size = 0.0f;
        row.final_proof_size = forward_proof_size;
        row.peak_memory = peak_mem;

        append_benchmark_csv("bench_results_infer.csv", row);
        std::fprintf(stderr, "[bench] appended to bench_results_infer.csv\n");

        return 0;
    }

    printf("RNN Model initialized.\n");

    vector<vector<vector<vector<F>>>> witness_matrixes_rnn(1);
    vector<vector<commitment>> comm_aggr_rnn(1);

    vector<struct proof> total_aggr_first_transcript;
    vector<F> total_aggr_first_hash;

    vector<vector<vector<vector<F>>>> proof_witness_matrixes_rnn(1);
    vector<vector<commitment>> proof_comm_aggr_rnn(1);

    float RNN_PoGDsize;

    float commitment_time = 0.0f;
    double recursion_verifier_time = 0.0;
    float aggregation_time_recursion = 0.0f;

    check_input_rnn(T,n);
    float rnn_data_prove_time = total_prove;
    printf("PoGD dataset proving time: %lf\n", rnn_data_prove_time);
    printf("PoGD dataset check size: %lf\n",proof_size(Transcript));

    F lr = quantize(0.01);

    Transcript.clear();

    vector<int> labels(T);
    vector<vector<F>> y_true(T, vector<F>(k, F(0)));
    for(int t = 0; t < T; ++t) {
        int cls = labels[t];        // 真正的类别编号
        y_true[t][cls] = quantize(1.0f);  // 量化域中的“1.0”
    }
	
	for(int i = 0; i < arity; i++){
		witness_matrixes_rnn[0].clear();
		comm_aggr_rnn[0].clear();

		vector<struct proof> total_pogd_transcript;
		vector<F> total_pogd_hash;
		for(int j = 0; j < arity; j++){
			printf("Iteration : %d\n", i*arity + j);

			rnn_net = rnn_forward(rnn_net,X_rnn,h0);
            DumpTrainingSnapshot(rnn_net, X_rnn, h0);
			rnn_net = rnn_bptt(rnn_net, y_true);
			rnn_net = rnn_update(rnn_net, lr);

			printf("RNN Compute finished.\n");

			printf("Begin model and witness commitment.\n");
			vector<F> witness_rnn;
			vector<F> new_model_rnn;

			get_witness_rnn(rnn_net, witness_rnn);
			printf("Witness size (before add model parameter): %zu\n", witness_rnn.size());
			get_model_rnn(rnn_net,new_model_rnn);
			printf("Model size : %zu\n", new_model_rnn.size());
			witness_rnn.insert(witness_rnn.end(),new_model_rnn.begin(),new_model_rnn.end());

			pad_vector(witness_rnn);
			printf("Total witness size (after padding):: %d\n", witness_rnn.size());
			
			vector<vector<F>> witness_matrix_rnn;
			commitment commitment_rnn;
			
			commitment_time = g_commit_general;
			poly_commit(witness_rnn, witness_matrix_rnn, commitment_rnn,levels);
			double model_commitment_time = g_commit_general - commitment_time;
			printf("Model Commit time : %lf\n", model_commitment_time);

			witness_matrixes_rnn[0].push_back(witness_matrix_rnn);
			witness_matrix_rnn.clear();
			comm_aggr_rnn[0].push_back(commitment_rnn);
			vector<vector<__hhash_digest>>().swap(commitment_rnn.hashes_sha);

			prove_rnn_update(rnn_net);
			prove_rnn_backward(rnn_net, h0);
			prove_rnn_forward(rnn_net,h0);

			adopt_next_and_drop_grads(rnn_net);
			RNN_PoGDsize = proof_size(Transcript);
			
			total_pogd_hash.insert(total_pogd_hash.end(), x_transcript.begin(), x_transcript.end());
			total_pogd_transcript.insert(total_pogd_transcript.end(), Transcript.begin(),Transcript.end());
			x_transcript.clear();
			Transcript.clear();
			double pogd_time = Forward_propagation_rnn + Backward_propagation_rnn+parameter_update_rnn;
			total_pogd += pogd_time;

			printf("PoGD proof size: %lf\n",RNN_PoGDsize);
			printf("Total PoGD GKR time : %lf\n",pogd_time-logup_acc_commit_time);
			printf("Total PoGD LogUp time : %lf\n", logup_acc_total_time);
			printf("Total PoGD Proving time : %lf\n",pogd_time);
		}

		printf("PoGD finished.\n");

		Transcript.clear();
		x_transcript.clear();


		double aggr_time_first;
		double aggr_now = total_aggregation_time;
		vector<F> proof_witness_rnn = prove_aggr(witness_matrixes_rnn, comm_aggr_rnn);
		printf("Aggregated Witness size : %zu\n", proof_witness_rnn.size());
		vector<proof> aggr_proofs_rnn = Transcript;
		vector<F> aggr_hashes_round = x_transcript;
		aggr_time_first = total_aggregation_time - aggr_now;
		aggregation_time_recursion += aggr_time_first;

		const double proof0   = total_prove;
		vector<struct proof> u_proof_round;
		aggr_proofs_rnn.insert(aggr_proofs_rnn.end(),total_pogd_transcript.begin(),total_pogd_transcript.end());
		//verify first round aggregation
		reset_gkr_verifier_time();
		u_proof_round.push_back(verify_proof(aggr_proofs_rnn));
		recursion_verifier_time += get_gkr_verifier_time();


		vector<proof> temp_P_rnn = mimc_sumcheck(total_pogd_hash);
		u_proof_round.insert(u_proof_round.end(), temp_P_rnn.begin(),temp_P_rnn.end());
		//第一轮aggregation中产生的hash的sumcheck
		temp_P_rnn = mimc_sumcheck(aggr_hashes_round);
		u_proof_round.insert(u_proof_round.end(), temp_P_rnn.begin(),temp_P_rnn.end());
		
		proof_witness_rnn.insert(proof_witness_rnn.end(),total_pogd_hash.begin(),total_pogd_hash.end());
		proof_witness_rnn.insert(proof_witness_rnn.end(), x_transcript.begin(),x_transcript.end());

		total_aggr_first_hash.insert(total_aggr_first_hash.end(), x_transcript.begin(),x_transcript.end());
		total_aggr_first_transcript.insert(total_aggr_first_transcript.end(), u_proof_round.begin(), u_proof_round.end());

		pad_vector(proof_witness_rnn);
		vector<vector<F>> proof_witness_matrix_rnn;
		commitment proof_com_rnn;
		printf("proof witness size : %d\n",proof_witness_rnn.size());
		const double proof1          = total_prove    - proof0;
		aggregation_time_recursion += proof1;
		
		commitment_time = g_commit_aggr;
		poly_commit(proof_witness_rnn, proof_witness_matrix_rnn, proof_com_rnn, levels, CommitSrc::Aggregation);
		double aggr_commitment_time = g_commit_aggr - commitment_time;
		printf("Aggregation Commit time : %lf\n", aggr_commitment_time);

		proof_witness_matrixes_rnn[0].push_back(proof_witness_matrix_rnn);
		proof_comm_aggr_rnn[0].push_back(proof_com_rnn);

		x_transcript.clear();
		Transcript.clear(); 
	}
	

	x_transcript.clear();
	Transcript.clear(); 

	double aggr_time_second;
	double aggr_now = total_aggregation_time;
	prove_aggr(proof_witness_matrixes_rnn,proof_comm_aggr_rnn);
	auto aggr_proofs_rnn_round2 = Transcript;
	vector<F> rnn_aggregation_second_hashes = x_transcript;
	x_transcript.clear();
	aggr_time_second = total_aggregation_time - aggr_now;
	aggregation_time_recursion += aggr_time_second;

	const double proof0  = total_prove;
	vector<proof> final_proofs;
	aggr_proofs_rnn_round2.insert(aggr_proofs_rnn_round2.end(), total_aggr_first_transcript.begin(),total_aggr_first_transcript.end());
	reset_gkr_verifier_time();
	final_proofs.push_back(verify_proof(aggr_proofs_rnn_round2));
	recursion_verifier_time += get_gkr_verifier_time();

	vector<proof> temp_P_rnn_second = mimc_sumcheck(total_aggr_first_hash);
	final_proofs.insert(final_proofs.end(),temp_P_rnn_second.begin(),temp_P_rnn_second.end());
	temp_P_rnn_second = mimc_sumcheck(rnn_aggregation_second_hashes);
	final_proofs.insert(final_proofs.end(),temp_P_rnn_second.begin(),temp_P_rnn_second.end());

	const double proof1          = total_prove     - proof0;
	aggregation_time_recursion += proof1;

	double wall0 = (double)clock();

	float final_proof_size = proof_size(final_proofs);

	printf("RNN_proof size: %lf\n",final_proof_size);

	reset_gkr_verifier_time();
	verify_proof(final_proofs);
	double final_verifier_time = get_gkr_verifier_time();
	double wall1 = (double)clock();
	aggregation_time_recursion += (wall1 - wall0) / CLOCKS_PER_SEC;

	printf("Average Update Prove time: %lf\n",total_update/((float)arity*(float)arity));printf("Average Backward Prove time: %lf\n",total_backward/((float)arity*(float)arity));
	printf("Average Forward Prove time: %lf\n",total_foward/((float)arity*(float)arity));
	printf("Average LogUp Prove time: %lf\n",total_logup_prove/((float)arity*(float)arity));
	printf("Average PoGD Prove time: %lf\n",total_pogd/((float)arity*(float)arity));


	printf("Total Recursion P : %lf, Amortized : %lf\n",aggregation_time_recursion,(aggregation_time_recursion)/((float)arity*(float)arity) );
	printf("Total Aggregation time : %lf, Amortized : %lf\n", total_aggregation_time,total_aggregation_time/((float)arity*(float)arity));
	printf("Total Aggregation proving time : %lf, Amortized : %lf\n", total_aggregation_prove_time,total_aggregation_prove_time/((float)arity*(float)arity));
	printf("Commitment time : %lf\n", (g_commit_general+g_commit_aggr+g_commit_logup)/((float)arity*(float)arity));
	printf("Verifier time: %lf\n", final_verifier_time);

	printf("Peak memory: %.2f GB\n", PeakMemGB());

	BenchmarkRow row;
	row.timestamp = now_iso8601();
	row.hostname = get_hostname();
	const char* rev    = std::getenv("GIT_REV");

	row.T = T; row.n = n; row.m = m; row.k = k;
	row.PC_scheme = PC_scheme;
	row.Commitment_hash = Commitment_hash;
    row.levels = levels; 
	row.arity = arity;
	row.param_count = 1LL*n*m + 1LL*m*m + 1LL*k*m + m + k;

	row.data_prove_time = rnn_data_prove_time;
	row.avg_update = total_update/((float)arity*(float)arity);
	row.avg_backward = total_backward/((double)arity*(double)arity);
    row.avg_forward = total_foward/((double)arity*(double)arity);
    row.avg_logup = total_logup/((double)arity*(double)arity);
	row.avg_logup_prove = total_logup_prove/((double)arity*(double)arity);
    row.avg_pogd = total_pogd/((double)arity*(double)arity);

	row.agg_recursion_total = aggregation_time_recursion;
    row.agg_recursion_amort = aggregation_time_recursion/((double)arity*(double)arity);
    row.aggr_time_total = total_aggregation_time;
	row.aggr_time_avg = total_aggregation_time/((float)arity*(float)arity);
    row.aggr_prove_total  = total_aggregation_prove_time;
	row.aggr_prove_avg = total_aggregation_prove_time/((float)arity*(float)arity);

	row.commit_time_avg = (g_commit_general+g_commit_aggr+g_commit_logup)/((float)arity*(float)arity);
	row.recursion_verifier_time = recursion_verifier_time/((float)arity*(float)arity);
	row.verifier_time = final_verifier_time;


	row.pogd_proof_size = RNN_PoGDsize;
	row.final_proof_size = final_proof_size;
	row.peak_memory = PeakMemGB();

	std::string out_path = "bench_results.csv";
	append_benchmark_csv(out_path, row);
    std::fprintf(stderr, "[bench] appended to %s\n", out_path.c_str());

	return 0;
}