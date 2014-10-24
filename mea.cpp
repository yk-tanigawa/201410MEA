#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

class sequence;
class mea;

typedef enum { M, X, Y, } tbltype_t;

class sequence{
  string header; 
  int length;
  vector<int> ary;
  int c_to_i(char);
  char i_to_c(int);
public:
  void set(string, string);
  void dump(int);
  int get_ary(int i){ return ary.at(i); }
  int len(){ return length; }
};
  
class mea{
  vector<sequence> seq;
  long double T, d, e;
  int alph_size;
  vector<vector <long double> > s;    /* substitution matrix */
  vector<vector<vector <long double> > > forward, backward;
  long double lpd, lpe;
  vector<vector <long double> > lpij;
  long double  lz;
  void read_data_seqs(ifstream &);
  void read_data_params(ifstream &);
  void forward_backward_init(int, int);
  void forward_calc(int, int);
  void backward_calc(int, int);
  void forward_backward_terminate(int, int);
public:
  void set_f(tbltype_t type, int i, int j, long double val){
    forward.at((int) type).at(i).at(j) = val;
  }
  void set_b(tbltype_t type, int i, int j, long double val){
    backward.at((int) type).at(i).at(j) = val;
  }
  long double get_f(tbltype_t type, int i, int j){
    return forward.at((int) type).at(i).at(j);
  }
  long double get_b(tbltype_t type, int i, int j){
    return backward.at((int) type).at(i).at(j);
  }
  long double get_lpij(int i, int j){
    return lpij.at(seq.at(0).get_ary(i)).at(seq.at(1).get_ary(j));
  }
  void read_data(char *, char *);
  void forward_backward();
  void dump(int);
};

template <class T>
void mat_dump(vector<vector <T> > matrix){
  for(int i = 0; i < matrix.size(); i++){
    for(int j = 0; j < matrix.at(i).size(); j++){
      cout << matrix.at(i).at(j) << " ";
    }
    cout << endl;
  }
  return;
}

istream &getline_wocomment(char c, istream &is, string &str){
  /*  read from stream without comments, * 
   *  'c' will be treated as a delimiter */
  getline(is, str);  int comment_start;
  while((comment_start= str.find(c)) == 0){getline(is, str);}
  if(comment_start > 0){ str.erase(comment_start); }
  return is;
}

inline void mea::forward_backward_init(int n, int m){
  /* forward variables と backward variablesのメモリを確保 */
  forward.resize(3); backward.resize(3);
  for(int t = 0; t < 3; t++){
    forward.at(t).resize(n);  backward.at(t).resize(n);
    for(int i = 0; i < n; i++){
      forward.at(t).at(i).resize(m, 0);
      backward.at(t).at(i).resize(m, 0);
    }
  }

  /* forward変数の初期化 */
  for(int i = 0; i < n; i++){
    set_f(M, i, 0, logl(0));
    set_f(X, i, 0, lpd + (i - 1) * lpe);
    set_f(Y, i, 0, logl(0)); 
  }
  for(int j = 0; j < m; j++){
    set_f(M, 0, j, logl(0));
    set_f(X, 0, j, logl(0));
    set_f(Y, 0, j, lpd + (j - 1) * lpe); 
  }
  set_f(M, 0, 0, logl(1)); 
  set_f(X, 0, 0, logl(0)); 
  set_f(Y, 0, 0, logl(0));

  /* backward変数の初期化 */  
  for(int i = 0; i < n - 1; i++){
    set_b(M, i, m - 1, lpd + (m - i - 1) * lpe);
    set_b(X, i, m - 1, (m - i) * lpe);
    set_b(Y, i, m - 1, logl(0)); 
  }
  for(int j = 0; j < m - 1; j++){
    set_b(M, n - 1, j, lpd + (n - j - 1) * lpe);
    set_b(X, n - 1, j, logl(0));
    set_b(Y, n - 1, j, (n - j) * lpe); 
  }
  set_b(M, n - 1, m - 1, logl(1)); 
  set_b(X, n - 1, m - 1, logl(0));
  set_b(Y, n - 1, m - 1, logl(0));

  return;
}

inline void mea::forward_calc(int n, int m){
  for(int i = 1; i < n; i++){
    for(int j = 1; j < m; j++){
      set_f(M, i, j, 
	    logl(expl( get_f(M, i - 1, j - 1) ) +
		 expl( get_f(X, i - 1, j - 1) ) +
		 expl( get_f(Y, i - 1, j - 1) )) +
	    get_lpij(i, j));
      set_f(X, i, j, 
	    logl(expl( get_f(M, i - 1, j) + lpd ) +
		 expl( get_f(X, i - 1, j) + lpe )));
      set_f(Y, i, j, 
	    logl(expl( get_f(M, i, j - 1) + lpd ) +
		 expl( get_f(Y, i, j - 1) + lpe )));
    }
  }
  return;
}

inline void mea::backward_calc(int n, int m){
  for(int i = n - 2; i >= 0; i--){
    for(int j = m - 2; j >= 0; j--){
      set_b(M, i, j, 
	    logl(expl(get_b(M, i + 1, j + 1) 
		      + get_lpij(i + 1, j + 1)) +
		 expl(get_b(X, i + 1, j) + lpd) +
		 expl(get_b(Y, i, j + 1) + lpd)));
      set_b(X, i, j, 
	    logl(expl(get_b(M, i + 1, j + 1) 
		      + get_lpij(i + 1, j + 1)) +
		 expl(get_b(X, i + 1, j) + lpe)));
      set_b(Y, i, j, 
	    logl(expl(get_b(M, i + 1, j + 1) 
		      + get_lpij(i + 1, j + 1)) +
		 expl(get_b(Y, i, j + 1) + lpe)));
    }
  }
  return;
}

inline void mea::forward_backward_terminate(int n, int m){
  lz =  logl(expl(get_f(M, n - 1, m - 1)) +
	     expl(get_f(X, n - 1, m - 1)) + 
	     expl(get_f(Y, n - 1, m - 1)));
  cout << lz << endl;
  cout << expl(lz) << endl;
  long double temp 
    = logl(expl(get_b(M, 0, 0) + get_lpij(0, 0)) + 
	   expl(get_b(X, 0, 0) + lpd) +
	   expl(get_b(Y, 0, 0) + lpd));
  cout << temp << endl;
  return;
}

inline void mea::forward_backward(){
  if( seq.size() < 2){ return; }
  int n = seq.at(0).len(), m = seq.at(1).len();
  forward_backward_init(n, m);
  forward_calc(n, m);
  backward_calc(n, m);
  forward_backward_terminate(n, m);
  return;
}


inline vector<string> split(const string &str, char delim){
  /* strを delim でsplitして vector<string> として返す */
  istringstream iss(str);  string tmp;  vector<string> res;
  while(getline(iss, tmp, delim)){ res.push_back(tmp); }
  return res;
}

int sequence::c_to_i(char c){
  switch(toupper(c)){
    case 'A': return 0; break;
    case 'C': return 1; break;
    case 'G': return 2; break;
    case 'T':
    case 'U': return 3; break;
    default: return -1;
  }
}

char sequence::i_to_c(int i){
  switch(i){
    case 0: return 'A'; break;
    case 1: return 'C'; break;
    case 2: return 'G'; break;
    case 3: return 'T'; break;
    default: return ' ';
  }
}

void sequence::set(string h, string s){
  header = h; length = s.length(); ary.resize(length, 0);
  for(int i = 0; i < length; i++){ ary.at(i) = c_to_i(s.at(i)); }
  return;
}

void mea::read_data_seqs(ifstream &fs){
  string buf, head, sequence_str; sequence new_seq;
  while(getline(fs, buf)){
    if(buf[0] == '>'){
      if(sequence_str.length() > 0){
	new_seq.set(head, sequence_str); seq.push_back(new_seq);
      }
      head = buf;
    }else{ sequence_str += buf; }
  }
  new_seq.set(head, sequence_str); seq.push_back(new_seq);
  return;
}

void mea::read_data_params(ifstream &fs){
  string buf;
  getline_wocomment('%', fs, buf);
  sscanf(buf.c_str(), "%Lf", &T);
  getline_wocomment('%', fs, buf);
  sscanf(buf.c_str(), "%Lf %Lf", &d, &e);
  getline_wocomment('%', fs, buf);
  sscanf(buf.c_str(), "%d", &alph_size);
  s.resize(alph_size);
  for(int i = 0; i < alph_size; i++){
    s.at(i).resize(alph_size, 0);
    getline_wocomment('%', fs, buf);
    vector<string> input_str = split(buf, ' ');
    for(int j = 0; j < alph_size; j++){
      sscanf(input_str.at(j).c_str(), "%Lf", &(s.at(i).at(j)));
    }
  }

  /* logをとったパラメータの初期化 */
  lpij.resize(alph_size);
  for(int i = 0; i < alph_size; i++){
    lpij.at(i).resize(alph_size, 0);
    for(int j = 0; j < alph_size; j++){
      lpij.at(i).at(j) = s.at(i).at(j) / T;
    }
  }
  lpd = -1.0 * d/T;  lpe = -1.0 * e/T;
  return;
}

void mea::read_data(char *param, char *data){
  ifstream param_fs(param);
  if(param_fs.fail())
    { cerr << "cannnot open params file" << endl; exit(1); }
  read_data_params(param_fs);
  param_fs.close();

  ifstream data_fs(data);
  if(data_fs.fail())
    { cerr << "cannot open data file" << endl; exit(1); }
  read_data_seqs(data_fs);
  data_fs.close();
  
  return;
}

void sequence::dump(int width = INT_MAX){
  /* width でsequenceを表示するときの，一行あたりの文字数を指定できる */
  cout << header << " (length = " << length << ")" << endl;
  int l;
  for(l = 0; l < (length + width - 1) / width - 1; l++){
    for(int i = l * width; i < (l + 1) * width; i++){ cout << i_to_c(ary.at(i)); }
    cout << endl;
  }
  for(int i = l * width; i < length; i++){ cout << i_to_c(ary.at(i)); }
  cout << endl; return;
}

void mea::dump(int width = INT_MAX){
  /* width でsequenceを表示するときの，一行あたりの文字数を指定できる */
  cout << "T = " << T << endl;
  cout << "d = " << d << ", "
       << "e = " << e << endl;
  cout << "size of alphabet = " << alph_size << endl;
  mat_dump(s);
  for(int i = 0; i < seq.size(); i++){ seq.at(i).dump(width); }
  return;
}

int main(){
  mea job;
  job.read_data((char *)"params.txt", (char *)"sequences.txt");
  job.dump(60);
  job.forward_backward();
  return 0;
}
