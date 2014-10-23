#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

class sequence;
class mea;

class sequence{
  string header; 
  int length;
  vector<int> ary;
  int c_to_i(char);
  char i_to_c(int);
public:
  void set(string, string);
  void dump(int);
};
  
class mea{
  vector<sequence> seq;
  long double T;
  long double d, e;
  int alph_size;
  vector<vector <long double> > s; /* substitution matrix */
  vector<vector<vector <long double> > > forward, backward;
  void read_data_seqs(ifstream &);
  void read_data_params(ifstream &);
public:
  void read_data(char *, char *);
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
  return 0;
}
