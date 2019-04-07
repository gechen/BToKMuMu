#ifndef Physics_Analysis_CutOptimize_H
#define Physics_Analysis_CutOptimize_H

#include <vector>
#include "TCut.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <math.h>

using namespace std;

class OptStyle{
 public:
  OptStyle(bool i_flag, TString i_cut, TString i_sign, Double_t i_ini, Double_t i_step, Int_t i_len);
  OptStyle(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step, Int_t i_len);
  OptStyle(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step);
  OptStyle(bool i_flag, TString i_cut, Double_t i_ini, Double_t i_step);
  OptStyle(){};
  void set_cut(TString i_cut){ str_cut = i_cut;}
  void set_sign(TString i_sign){ str_sign = i_sign;}
  void set_ini(Double_t i_ini){ m_ini = i_ini;}
  void set_step(Double_t i_step){ m_step = i_step;}
  void set_len(Int_t i_len){ m_len = i_len;}
  void set_flag(bool i_flag){ m_flag = i_flag;}
  void cal_cut(Double_t i_via);
  void operator--(){ m_ini -= m_step; cal_cut(0);}
  void operator++(){ m_ini += m_step; cal_cut(0);}
  TCut get_cut(){ return m_cut;}
  bool get_flag(){ return m_flag;}
  void print(){ cout <<"mid:" << m_mid << " ini:" << m_ini << " step:" << m_step << "  cut:" <<  m_cut.GetTitle() << endl; }
 private:
  bool m_flag;
  TString str_abs, str_cut, str_sign;
  Double_t m_mid, m_ini, m_step;
  Int_t m_len;
  TCut m_cut;
  Int_t get_len(){ return m_len;}
};

OptStyle::OptStyle(bool i_flag, TString i_cut, Double_t i_ini, Double_t i_step):
  m_flag(i_flag),
     str_abs(""),
     str_cut(i_cut),
     str_sign(""),
     m_mid(0),
     m_ini(i_ini),
     m_step(i_step),
     m_len(1){
  cal_cut(0);
}


OptStyle::OptStyle(bool i_flag, TString i_cut, TString i_sign, Double_t i_ini, Double_t i_step, Int_t i_len):
  m_flag(i_flag),
     str_abs(""),
     str_cut(i_cut),
     str_sign(i_sign),
     m_mid(0),
     m_ini(i_ini),
     m_step(i_step),
     m_len(i_len){
  cal_cut(0);
}

OptStyle::OptStyle(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step, Int_t i_len):
   m_flag(i_flag),
     str_abs(i_fabs),
     str_cut(i_cut),
     str_sign(i_sign),
     m_mid(i_mid),
     m_ini(i_ini),
     m_step(i_step),
     m_len(i_len){
  cal_cut(0);
}


OptStyle::OptStyle(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step):
     m_flag(i_flag),
     str_abs(i_fabs),
     str_cut(i_cut),
     str_sign(i_sign),
     m_mid(i_mid),
     m_ini(i_ini),
     m_step(i_step),
     m_len(1){
  cal_cut(0);
}

void OptStyle::cal_cut(Double_t i_via){
  Double_t cur_value = m_ini + i_via*m_step;
  TString temp_str_cut("");
  if(m_len>1){
    if(str_abs==""){
      for(Int_t i=0; i<m_len; i++){
	temp_str_cut += str_cut;
	temp_str_cut += "[";
	temp_str_cut += i;
	temp_str_cut += "]";
	temp_str_cut += str_sign;
	temp_str_cut += cur_value;
	if(i<m_len-1) temp_str_cut += " && ";
      }
    } // end of abs if
    else {
      for(Int_t i=0; i<m_len; i++){
	temp_str_cut += str_abs;
	temp_str_cut += "(";
	temp_str_cut += str_cut;
	temp_str_cut += "[";
	temp_str_cut += i;
	temp_str_cut += "] - ";
	temp_str_cut += m_mid;
	temp_str_cut += ")";
	temp_str_cut += str_sign;
	temp_str_cut += cur_value;
	if(i<m_len-1) temp_str_cut += " && ";
      }
    } // end of abs if else 
  } // end of len > 1
  else{
    if(str_abs==""){
      temp_str_cut = str_cut + str_sign;
      temp_str_cut += cur_value;
    }
    else{
      temp_str_cut = str_abs;
      temp_str_cut += "(";
      temp_str_cut += str_cut;
      temp_str_cut += "-";
      temp_str_cut += m_mid;
      temp_str_cut +=  ")";
      temp_str_cut += str_sign;
      temp_str_cut += cur_value;
    } // end of abs if else 
  } // end of len > 1 if else 
  m_cut = temp_str_cut;
}


//===============================
//*****************************
//===============================
class CutOptimize{
 public:
  CutOptimize();
  void add_tree_sig(TTree *i_tree_sig, Double_t i_scale=1.0){ tree_sig.push_back(i_tree_sig); scale_sig.push_back(i_scale);}
  void set_tree_sig(TTree *i_tree_sig, Double_t i_scale=1.0){ tree_sig.clear(); scale_sig.clear(); add_tree_sig(i_tree_sig, i_scale); }
  void add_tree_sigbg(TTree *i_tree_sigbg, Double_t i_scale=1.0){ tree_sigbg.push_back(i_tree_sigbg); scale_sigbg.push_back(i_scale);}
  void set_tree_sigbg(TTree *i_tree_sigbg, Double_t i_scale=1.0){ tree_sigbg.clear(); scale_sigbg.clear(); add_tree_sigbg(i_tree_sigbg, i_scale); }
  void add_opt(OptStyle &i_sty){ v_opt.push_back(i_sty);}
  void add_opt(bool i_flag, TString i_cut, Double_t i_ini, Double_t i_step);
  void add_opt(bool i_flag, TString i_cut, TString i_sign, Double_t i_ini, Double_t i_step,Int_t i_len);
  void add_opt(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step, Int_t i_len);
  void add_opt(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step);
  Double_t cal_fom();
  Double_t cal_fom(Int_t i, Double_t f); // i: index in loop; f: sign of step
  void set_debug(bool i_debug){ m_debug = i_debug; }
  void loop();
  void print();
  void set_hard_cut(const char *i_cut){ m_hard_cut = i_cut; }
  void set_hard_cut(TCut i_cut){ m_hard_cut = i_cut; }
 private:
  vector<TTree *> tree_sig, tree_sigbg;
  vector<Double_t> scale_sig, scale_sigbg;
  vector<OptStyle> v_opt;
  TCut m_cur_cut, m_hard_cut;
  Double_t m_old_fom, m_new_fom;
  void print(Int_t i);
  bool m_debug;
};

CutOptimize::CutOptimize():
  tree_sig(0),
     tree_sigbg(0),
     m_old_fom(0),
     m_new_fom(0){
  v_opt.clear();
}

void CutOptimize::add_opt(bool i_flag, TString i_cut, Double_t i_ini, Double_t i_step){
  OptStyle temp_opt(i_flag, i_cut, i_ini, i_step);
  v_opt.push_back(temp_opt);
}

void CutOptimize::add_opt(bool i_flag, TString i_cut, TString i_sign, Double_t i_ini, Double_t i_step, Int_t i_len)
{
  OptStyle temp_opt(i_flag, i_cut, i_sign, i_ini, i_step, i_len);
  v_opt.push_back(temp_opt);
}

void CutOptimize::add_opt(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step, Int_t i_len){
  OptStyle temp_opt(i_flag, i_fabs, i_cut, i_sign, i_mid, i_ini, i_step, i_len);
  v_opt.push_back(temp_opt);
}

void CutOptimize::add_opt(bool i_flag, TString i_fabs, TString i_cut, TString i_sign, Double_t i_mid, Double_t i_ini, Double_t i_step){
  OptStyle temp_opt(i_flag, i_fabs, i_cut, i_sign, i_mid, i_ini, i_step);
  v_opt.push_back(temp_opt);
}


Double_t CutOptimize::cal_fom(){
  if(m_debug) cout << "in cal_fom()" << endl;
  m_cur_cut = m_hard_cut;
  for(Int_t i=0; i<(Int_t)v_opt.size(); i++) m_cur_cut += v_opt[i].get_cut();
  if(m_debug) cout << m_cur_cut.GetTitle() << endl;
  Double_t temp_sig_num = 0;
  for(int j=0; j<(Int_t)tree_sig.size(); j++){
    temp_sig_num += (tree_sig[j]->GetEntries(m_cur_cut))*(scale_sig[j]);
    if(m_debug) cout << "index of sig: " << j << "  sig num: " << temp_sig_num << endl;
  }
  Double_t temp_sigbg_num = 0; 
  for(int k=0; k<(Int_t)tree_sigbg.size(); k++){
    if(m_debug) cout << "tree address: " << tree_sigbg[k] << endl;
    temp_sigbg_num += (tree_sigbg[k]->GetEntries(m_cur_cut))*(scale_sigbg[k]);
    if(m_debug) cout << "index of sigbg: " << k << "  sigbg num: " << temp_sigbg_num << endl;
  }
  //  if(m_debug) cout << "inn cal_fom() before return, the fom: " << temp_sig_num/sqrt(temp_sigbg_num) << endl;
  if(temp_sigbg_num) return temp_sig_num/sqrt(temp_sigbg_num);
  else return 0;
}

Double_t CutOptimize::cal_fom(Int_t i_i, Double_t i_f){
  TCut temp_cur_cut = m_hard_cut;
  for(Int_t i=0; i<(Int_t)v_opt.size(); i++){
    if(i==i_i){
      v_opt[i].cal_cut(i_f);
      temp_cur_cut += v_opt[i].get_cut();
    }
      else temp_cur_cut += v_opt[i].get_cut();
  }
  Double_t temp_sig_num = 0;
  for(int j=0; j<(Int_t)tree_sig.size(); j++) temp_sig_num += (tree_sig[j]->GetEntries(temp_cur_cut))*(scale_sig[j]);
  Double_t temp_sigbg_num = 0; 
  for(int k=0; k<(Int_t)tree_sigbg.size(); k++) temp_sigbg_num += (tree_sigbg[k]->GetEntries(temp_cur_cut))*(scale_sigbg[k]);
  if(temp_sigbg_num) return temp_sig_num/sqrt(temp_sigbg_num);
  else return 0;
}

void CutOptimize::loop(){
  m_new_fom = cal_fom();
  while(m_new_fom>m_old_fom){
    m_old_fom = m_new_fom;
    for(Int_t i=0; i<(Int_t)v_opt.size(); i++){
      if(!(v_opt[i].get_flag())){ print(i); continue; }
      Double_t m_temp_fom(0);
      while(m_new_fom>m_temp_fom){
	m_temp_fom = m_new_fom;
	Double_t temp_low_fom = cal_fom(i,-1);
	Double_t temp_high_fom = cal_fom(i,1);
	Double_t temp_ori_fom = cal_fom(i,0);
	//if(m_debug) cout << "debug, low fom= " << temp_low_fom << "  ori fom= " << temp_ori_fom << "  high fom= " << temp_high_fom << endl;
	if(temp_low_fom > temp_ori_fom && temp_low_fom > temp_high_fom){
	  --(v_opt[i]);
	  m_new_fom = temp_low_fom;
	}
	if(temp_high_fom > temp_ori_fom && temp_high_fom > temp_low_fom){
	  ++(v_opt[i]);
	  m_new_fom = temp_high_fom;
	}
      } // end of one parameter optimization
      print(i);
    } // end of once vector loop
  }// end of whole loop
}

void CutOptimize::print(){
  cout << "FOM: " << m_new_fom << endl;
  cout << m_cur_cut.GetTitle() << endl << endl;
}

void CutOptimize::print(Int_t i_i){
  cal_fom();
  cout << "FOM: " << m_new_fom << endl;
  v_opt[i_i].print();
}

#endif
