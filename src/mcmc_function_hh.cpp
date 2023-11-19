#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <RcppParallel.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp ;
using namespace RcppParallel;

//########################################################################################################################################
//function for computing the serial interval density
// [[Rcpp::export]]
NumericVector serial_density(double p1,
double p2){
NumericVector a=NumericVector::create(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0);
return (pweibull(a+1,p1,p2)-pweibull(a,p1,p2));
}

//########################################################################################################################################
// function to generate normal random variable

double rnorm(double a, double b) { // a: mean, b: s.d.
	double c=a+b*sum(rnorm(1));
    return c;
}

//########################################################################################################################################
// function to generate binomial random number

int gen_binom(double p){
double cut=(double)rand()/(RAND_MAX);
int out=0;
if (cut<p){
out=1;
}
return out;
}

//########################################################################################################################################
// function to general multiviariable normal given sigma

NumericVector rmnorm(arma::mat sigma) {
int ncols=sigma.n_cols;
arma::rowvec c=arma::randn(1,ncols);
arma::rowvec a=c*arma::chol(sigma);   
NumericVector b=NumericVector(a.begin(),a.end());   
return b;
}

//########################################################################################################################################
//function to compute the prior likelihood 
// [[Rcpp::export]]
double prior_loglik(NumericVector para){
// check if the para are within their possible range
NumericVector out(para.length());
int b1;

//out(0)=R::dgamma(pow(1.0/para(0),2.0),1.5,1/0.0001,1); 
out(0)=R::dunif(para(0),0.009,5.00,1);
for (b1=2;b1>=1;--b1){
out(b1)=R::dunif(para(b1),0.000000000000000001,9.99,1);
}

out(3)=R::dunif(para(3),0.0,1.0,1);


for (b1=para.length()-1;b1>=4;--b1){
out(b1)=R::dnorm(para(b1),0.0,3.0,1);
}

double output=sum(out);
// if the prior is outside the parameter space
if (output< -9999999){
output=-9999999;
}
return output;
}




//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel for doing simulation
struct SimData:public Worker{
// source vector
RMatrix<double> data11;
RMatrix<double> data1;
RVector<double> para;
RVector<double> SI;
int n_inf;
int n_sus;
int with_rm;
int sep1;
int sep2;
RMatrix<int> record;
// destination vector
// initialize with source and destination
SimData(NumericMatrix data11,
NumericMatrix data1,
NumericVector para,
NumericVector SI,
int n_inf,
int n_sus,
int with_rm,
int sep1,
int sep2,
IntegerMatrix record) 
:data11(data11),data1(data1),para(para),SI(SI),n_inf(n_inf),n_sus(n_sus),with_rm(with_rm),sep1(sep1),sep2(sep2),record(record){}
void operator()(std::size_t begin, std::size_t end) {

// section to write the parallel version
// functor (pass input and output matrixes)

for (unsigned int b1=begin;b1<end;++b1){	
int b2;
int b3;
int b4;
int b5;
int b6;
double hazard;
double sus[int(data11(b1,1))];

// first set every household contacts have 0 infection status
for (b2=data11(b1,1)-1;b2>=1;--b2){
// initial point to set all people is uninfected at the beginning	
data11(b1,b2*sep2+sep1)=0;
data11(b1,b2*sep2+sep1+1)=-1;

// data
//## 1. hhID, 2.size, 3 influenza type ,4-5: follow start and end date 
//## 5 + 1. infectious status, 2. infection time, 3,4.start and end of follow-up 5.id 
//## 6.age, 7.sex, 8.vaccination, 9.oseltamivir treatment, 10.community type
sus[b2]=0;
if (n_sus>0){
for (b4=n_sus-1;b4>=0;--b4){
sus[b2]+=para[4+n_inf+b4]*(data11(b1,b3*sep2+sep1+3+n_inf+b4));	
}
}
}


for (b2=data11(b1,1)-1;b2>=0;--b2){
if (with_rm==1){
// generate the infectivity 
data11(b1,b2*sep2+sep1+2)=R::rnorm(0.0,para[0]);	
}
else{
data11(b1,b2*sep2+sep1+2)=0;	
}
}

// first need to have a time index
for (b2=data11(b1,3)+1;b2<=data11(b1,4);++b2){
// participant index
// index case do not need update
for (b3=data11(b1,1)-1;b3>=1;--b3){
if (data11(b1,b3*sep2+sep1)==0){
// compute the community hazard
hazard=para[1];
// compute the household hazard
for (b4=data11(b1,1)-1;b4>=0;--b4){
if ((b4!=b3)&&(data11(b1,b4*sep2+sep1)==1)){
// need to with the range of serial interval to have contribution
if ((b2-data11(b1,b4*sep2+sep1+1)>0)&&(b2-data11(b1,b4*sep2+sep1+1)<=SI.length())){
double hrisk=para[2];	
double inf=0;
// here need to add factor affecting infectivity
if (n_inf>0){
for (b6=n_inf-1;b6>=0;--b6){
inf+=para[4+b6]*(data11(b1,b4*sep2+sep1+3+b6));	
}
}
hrisk*=exp(inf);
hazard+=hrisk*SI[b2-data11(b1,b4*sep2+sep1+1)-1]*exp(data11(b1,b4*sep2+sep1+2)*(b4>=0))/pow(data11(b1,1)-1.0,para[3]);; // need to -1 beacuse the index is from 0 to 9
}
}	
}
// then need to add factor affecting transmission
hazard*=exp(sus[b3]);
// generate infeciton
// infection
if (gen_binom(1-exp(-hazard))){	
data11(b1,b3*sep2+sep1)=1;
data11(b1,b3*sep2+sep1+1)=b2;
}
}
}
}




}
}
};

//########################################################################################################################################
//function to do simulation
// [[Rcpp::export]]
List sim_data(NumericMatrix data1,
NumericVector SI,
NumericVector para,
int n_inf,
int n_sus,
int with_rm,
int sep1,      // sep1=5
int sep2){     // sep2=10
int b1;
int b2;
int b3;
// clone the data first
NumericMatrix data11(clone(data1));

IntegerMatrix record(data11.nrow(),100);


//1. hhID, 2.size, 3 influenza type ,4-5: follow start and end date 
//5 + 1. infectious status, 2. infection time, 3,4.start and end of follow-up 5.id 
//6.age, 7.sex, 8.vaccination, 9.oseltamivir treatment, 10.community type

// compute the serial_density to use
//NumericVector SI=serial_density(para[0],para[1]);

// call parallel program
SimData simdata1(data11,data1,para,SI,n_inf,n_sus,with_rm,sep1,sep2,record);

// call parallelFor to do the work
parallelFor(0,data1.nrow(),simdata1);


return List::create(_[""]=data11,
_[""]=record);
} 


//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel for Digraphlikelihood
struct LogLik:public Worker{
// source vector
RMatrix<double> out;
RMatrix<double> out2;
RMatrix<double> data;
RVector<double> para;
RVector<double> SI;
int n_inf;
int n_sus;
int with_rm;
int sep1;
int sep2;
RMatrix<double> record;
// destination vector
// initialize with source and destination
LogLik(NumericMatrix out,
NumericMatrix out2,	
NumericMatrix data,
NumericVector para,
NumericVector SI,
int n_inf,
int n_sus,
int with_rm,
int sep1,
int sep2,
NumericMatrix record) 
:out(out),out2(out2),data(data),para(para),SI(SI),n_inf(n_inf),n_sus(n_sus),with_rm(with_rm),sep1(sep1),sep2(sep2),record(record){}
void operator()(std::size_t begin, std::size_t end) {

// section to write the parallel version
// functor (pass input and output matrixes)
for (unsigned int b1=begin;b1<end;++b1){
int b2;
int b3;
int b4;
int b5;
int b6;

// b3 is the participant index
for (b3=data(b1,1)-1;b3>=1;--b3){
//if (!((data(b1,b3*sep2+sep1)==1)&(data(b1,b3*sep2+sep1+1)==data(b1,sep1+1)))){
// need to add factor addecting susceptibility
double sus=0;
if (n_sus>0){
for (b4=n_sus-1;b4>=0;--b4){
sus+=para[4+n_inf+b4]*(data(b1,b3*sep2+sep1+3+n_inf+b4));	
}
}
//sus+=para[3]*(data(b1,1)==2);
//sus+=para[6]*(data(b1,b3*sep2+sep1+5)<=5);


// the final date for contribution from non-infection
int finaltime=data(b1,4)+1;
if (data(b1,b3*sep2+sep1)==1){
finaltime=data(b1,b3*sep2+sep1+1);	
}

double h[finaltime-int(data(b1,3))];
// fill the community risk
for (b2=finaltime-data(b1,3)-1;b2>=0;--b2){
h[b2]=para[1];	
}


for (b4=data(b1,1)-1;b4>=0;--b4){
if ((b4!=b3)&&(data(b1,b4*sep2+sep1)==1)){
for (b5=SI.length()-1;b5>=0;--b5){
//if (data(b1,b4*sep2+sep1+1)+b5<=data(b1,4)){
// if infection date of individual b3 is smaller than the final time	
if (data(b1,b4*sep2+sep1+1)+b5+1<=finaltime){ 
double hrisk=para[2];
double inf=0;
// here need to add factor affecting infectivity
if (n_inf>0){
for (b6=n_inf-1;b6>=0;--b6){
inf+=para[4+b6]*(data(b1,b4*sep2+sep1+3+b6));	
}
}
hrisk*=exp(inf);
h[int(data(b1,b4*sep2+sep1+1)-data(b1,3))+b5]+=hrisk*SI[b5]*exp(data(b1,b4*sep2+sep1+2)*(b4>=0))/pow(data(b1,1)-1.0,para[3]);;
}
}
}	
}
//}


for (b2=finaltime-data(b1,3)-2;b2>=0;--b2){
out(b1,b3)-=h[b2]*exp(sus);	
}
if (data(b1,b3*sep2+sep1)==1){
out(b1,b3)+=log(1-exp(-h[finaltime-int(data(b1,3))-1]*exp(sus)));	
}

//if (b3==1){
//for (b2=finaltime-data(b1,3)-1;b2>=0;--b2){
//record(b1,b2)=h[b2];
//}
//}

//}


}

if (with_rm){
for (b3=data(b1,1)-1;b3>=0;--b3){
// input the rm for infectivity
if (data(b1,b3*sep2+sep1)==1){
out2(b1,b3)=R::dnorm(data(b1,b3*sep2+sep1+2),0.0,para[0],1);
}
}
}

}
}
};

//########################################################################################################################################
//function to likelihood
// [[Rcpp::export]]
List loglik(NumericMatrix data,
NumericVector SI,
NumericVector para,
int n_inf,
int n_sus,
int with_rm,
int sep1,
int sep2){
int max_member=max(data(_,1));
// check if the para are within their possible range
NumericMatrix out(data.nrow(),max_member);
NumericMatrix out2(data.nrow(),max_member);
NumericMatrix record(1,1);
//NumericVector SI=serial_density(para[0],para[1]);
// call parallel program
LogLik loglik(out,out2,data,para,SI,n_inf,n_sus,with_rm,sep1,sep2,record);
// call parallelFor to do the work
parallelFor(0,data.nrow(),loglik);


return List::create(_[""]=out,
	_[""]=out2,
	_[""]=record);

}




//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel function
struct AllUpdate:public Worker{
// source vector
RMatrix<double> dataout;
RMatrix<double> datapro;
RMatrix<double> dataorg;
RMatrix<double> loglik1out;
RMatrix<double> loglik2out;
RMatrix<double> loglik1pro;
RMatrix<double> loglik2pro;
RMatrix<double> mcmcrecord;
RMatrix<double> data;
RVector<double> SI;
RVector<double> para;
int member;
RMatrix<double> loglik1;
RMatrix<double> loglik2;
RMatrix<double> temprecord;
int n_inf;
int n_sus;
int with_rm;
int sep1;
int sep2;
// destination vector
// initialize with source and destination
AllUpdate(NumericMatrix dataout,
NumericMatrix datapro,
NumericMatrix dataorg,
NumericMatrix loglik1out,
NumericMatrix loglik2out,
NumericMatrix loglik1pro,
NumericMatrix loglik2pro,
NumericMatrix mcmcrecord,
NumericMatrix data,
NumericVector SI,
NumericVector para,
int member,
NumericMatrix loglik1,
NumericMatrix loglik2,
NumericMatrix temprecord,
int n_inf,
int n_sus,
int with_rm,
int sep1,
int sep2) 
:dataout(dataout),datapro(datapro),dataorg(dataorg),loglik1out(loglik1out),loglik2out(loglik2out),loglik1pro(loglik1pro),loglik2pro(loglik2pro),mcmcrecord(mcmcrecord),data(data),SI(SI),para(para),member(member),loglik1(loglik1),loglik2(loglik2),temprecord(temprecord),n_inf(n_inf),n_sus(n_sus),with_rm(with_rm),sep1(sep1),sep2(sep2){}

void operator()(std::size_t begin, std::size_t end) {

// section to write the parallel version
// functor (pass input and output matrixes)
for (unsigned int b1=begin;b1<end;++b1){
int b2;
int b3;
int b4;
int b5;
int b6;

if (member<data(b1,1)){
if ((with_rm==1)|((member!=0)&(datapro(b1,member*sep2+sep1+1)!=-1))){
if (data(b1,member*sep2+sep1)==1){

// first update infection time

double proratio=0;


// the get the propose infection time
// everyone needs to update infection time

if (member!=0){
if (dataorg(b1,member*sep2+sep1+1)==-1){
datapro(b1,member*sep2+sep1+1)=rand()%int((data(b1,4)-data(b1,3)))+data(b1,3);
}
}

//if (member==0){
if (with_rm==1){
datapro(b1,member*sep2+sep1+2)=data(b1,member*sep2+sep1+2)+R::rnorm(0.0,para[0]);
}
//}
//proratio+=log(ILI(data1pro(b1,5+5*member+1)-1,data1(b1,46)))-log(ILI(data1(b1,5+5*member+1)-1,data1(b1,46)));





// here compute the likelihood for propose
// first level is transmission

// b3 is the participant index
for (b3=datapro(b1,1)-1;b3>=1;--b3){
loglik1pro(b1,b3)=0;	
//if (!((datapro(b1,b3*sep2+sep1)==1)&(datapro(b1,b3*sep2+sep1+1)==datapro(b1,sep1+1)))){
// need to add factor addecting susceptibility
double sus=0;
//sus+=para[3]*(datapro(b1,1)==2);
//sus+=para[6]*(datapro(b1,b3*sep2+sep1+5)<=5);
if (n_sus>0){
for (b4=n_sus-1;b4>=0;--b4){
sus+=para[4+n_inf+b4]*(datapro(b1,b3*sep2+sep1+3+n_inf+b4));	
}
}

// the final date for contribution from non-infection
int finaltime=datapro(b1,4)+1;
if (datapro(b1,b3*sep2+sep1)==1){
finaltime=datapro(b1,b3*sep2+sep1+1);	
}

double h[finaltime-int(datapro(b1,3))];
// fill the community risk
for (b2=finaltime-datapro(b1,3)-1;b2>=0;--b2){
h[b2]=para[1];	
}


for (b4=datapro(b1,1)-1;b4>=0;--b4){
if ((b4!=b3)&&(datapro(b1,b4*sep2+sep1)==1)){
for (b5=SI.length()-1;b5>=0;--b5){
//if (datapro(b1,b4*sep2+sep1+1)+b5<=datapro(b1,4)){
// if infection date of individual b3 is smaller than the final time	
if (datapro(b1,b4*sep2+sep1+1)+b5+1<=finaltime){ 
double hrisk=para[2];
double inf=0;
// here need to add factor affecting infectivity
if (n_inf>0){
for (b6=n_inf-1;b6>=0;--b6){
inf+=para[4+b6]*(datapro(b1,b4*sep2+sep1+3+b6));	
}
}
hrisk*=exp(inf);
h[int(datapro(b1,b4*sep2+sep1+1)-datapro(b1,3))+b5]+=hrisk*SI[b5]*exp(datapro(b1,b4*sep2+sep1+2)*(b4>=0))/pow(datapro(b1,1)-1.0,para[3]);;
}
}
}	
}
//}


for (b2=finaltime-datapro(b1,3)-2;b2>=0;--b2){
loglik1pro(b1,b3)-=h[b2]*exp(sus);	
}
if (datapro(b1,b3*sep2+sep1)==1){
loglik1pro(b1,b3)+=log(1-exp(-h[finaltime-int(datapro(b1,3))-1]*exp(sus)));	
}

}

// second level is the infectivity random effect
if (with_rm==1){
loglik2pro(b1,member)=R::dnorm(datapro(b1,member*sep2+sep1+2),0.0,para[0],1);
}



// here do the metropolis hasting update
double liknew=loglik2pro(b1,member);
double likold=loglik2(b1,member);

// don't add level2, beacuse using gibbs sampler
for (b2=loglik1pro.ncol()-1;b2>=0;--b2){
liknew+=loglik1pro(b1,b2);
likold+=loglik1(b1,b2);	
}


double loglikratio=liknew-likold;
double accept_pro=pow(exp(1),loglikratio-proratio);
mcmcrecord(b1,1)=accept_pro;
mcmcrecord(b1,2)=loglikratio;	
mcmcrecord(b1,3)=proratio;
if (gen_binom(accept_pro)){
mcmcrecord(b1,0)=1;	
// if accept, make the out to the same as the proposal
dataout(b1,member*sep2+sep1+1)=datapro(b1,member*sep2+sep1+1);
dataout(b1,member*sep2+sep1+2)=datapro(b1,member*sep2+sep1+2);
for (b2=loglik1pro.ncol()-1;b2>=0;--b2){
loglik1out(b1,b2)=loglik1pro(b1,b2);	
}
// add the lik for rm
loglik2out(b1,member)=loglik2pro(b1,member);
}
else{
mcmcrecord(b1,0)=-1;	
}

}

}
}

}
}
};




//########################################################################################################################################
// function to update infection time, baseline AT titer
// but not the infection status
// [[Rcpp::export]]
List all_update(NumericMatrix data,
NumericMatrix dataorg,	
NumericVector SI,
NumericVector para,
int member,
NumericMatrix loglik1,
NumericMatrix loglik2,
int n_inf,
int n_sus,
int with_rm,
int sep1,
int sep2){

NumericMatrix datapro(clone(data));
NumericMatrix dataout(clone(data));
NumericMatrix loglik1pro(clone(loglik1));
NumericMatrix loglik2pro(clone(loglik2));
NumericMatrix loglik1out(clone(loglik1));
NumericMatrix loglik2out(clone(loglik2));
// 1.accept/reject
NumericMatrix mcmcrecord(data.nrow(),20);
NumericMatrix temprecord(data.nrow(),20);

int b1;


// call parallel program
AllUpdate allupdate(dataout,datapro,dataorg,loglik1out,loglik2out,loglik1pro,loglik2pro,mcmcrecord,data,SI,para,member,loglik1,loglik2,temprecord,n_inf,n_sus,with_rm,sep1,sep2);
// call parallelFor to do the work
parallelFor(0,data.nrow(),allupdate);


double a1=0;
double a2=0;
for (b1=mcmcrecord.nrow()-1;b1>=0;--b1){
if (mcmcrecord(b1,0)==1){
++a1;
}
if (mcmcrecord(b1,0)!=0){
++a2;	
}
}


return List::create(_[""]=dataout,
_[""]=loglik1out,
_[""]=loglik2out,
_[""]=a1/a2,
_[""]=mcmcrecord,
_[""]=datapro,
_[""]=loglik1pro,
_[""]=loglik2pro,
_[""]=temprecord);
}




//##############################################################################################################################################
//##############################################################################################################################################
// function for mcmc
// [[Rcpp::export]]
List mcmc(NumericMatrix data1,
NumericVector SI,
int mcmc_n,             // length of mcmc stain
int burnin,
int thinning,
NumericVector int_para, // initial parameter
NumericVector move,     // which one should move in the model
NumericVector sigma,
int n_inf,
int n_sus,
int with_rm,
int sep1,
int sep2){            

// create the vector for use
int b0;
int b1;
int b2;
int b3;
int b4;
int moveindex;
int max_member=max(data1(_,1));
//####################################################################################################################################
// backup data 
NumericMatrix data11(clone(data1));



// matrix to record LL
// need to set number of parameter here
NumericMatrix p_para(mcmc_n,int_para.length());
NumericMatrix p_para_r(mcmc_n,sum(move));
p_para(0,_)=int_para;
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(0,moveindex)=p_para(0,b1);
--moveindex;
}	
}

// first row is the overall matrix, other three row is the indiviudal likelihood
NumericVector acceptrate(int_para.length());
NumericMatrix LL1(mcmc_n,3);
NumericMatrix LL2(mcmc_n,3);

//####################################################################################################################################
// initial step

//####################################################################################################################################
// compute likelihood
List loglikall=loglik(data11,SI,p_para(0,_),n_inf,n_sus,with_rm,sep1,sep2);
List loglikallpro;
NumericMatrix loglik1=loglikall(0);
NumericMatrix loglik2=loglikall(1);
NumericMatrix loglik1pro;
NumericMatrix loglik2pro;

LL1(0,1)=sum(loglik1);
LL1(0,2)=sum(loglik2);
LL1(0,0)=LL1(0,1)+LL1(0,2);

NumericVector temploglik(3);
NumericVector newloglik(3);
temploglik(0)=LL1(0,0)+prior_loglik(p_para(0,_));
temploglik(1)=LL1(0,1);
temploglik(2)=LL1(0,2);


double loglikeratio;
double accept_pro;
NumericVector pro_para(int_para.length());



//####################################################################################################################################
// main mcmc step
NumericMatrix updateacceptrate(mcmc_n,20);

int rownumber=floor((mcmc_n-burnin)/thinning);
// record the random effect
NumericMatrix rmrecord(rownumber,data11.nrow());
int recindex=0;
//####################################################################################################################################
for (b0=1;b0<mcmc_n;++b0){

// after 500 step, then set the sigma to be the empirical sigma
if (b0>500){
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){	
NumericVector temp1(b0-1);
for (b2=b0-2;b2>=0;--b2){
temp1(b2)=p_para(b2,b1);	
}
sigma(b1)=sd(temp1);
// tuning
if (acceptrate(b1)<0.1){
sigma(b1)*=0.5;
}	
if ((acceptrate(b1)<0.15)&(acceptrate(b1)>0.1)){
sigma(b1)*=0.8;
}
if ((acceptrate(b1)<0.2)&(acceptrate(b1)>0.15)){
sigma(b1)*=0.95;
}
if ((acceptrate(b1)<0.4)&(acceptrate(b1)>0.3)){
sigma(b1)*=1.05;
}
if ((acceptrate(b1)<0.9)&(acceptrate(b1)>0.4)){
sigma(b1)*=1.2;
}
if (acceptrate(b1)>0.9){
sigma(b1)*=2;
}
}
}
}

// metorpolis-hasing update on parameter
for (b1=0;b1<int_para.length();++b1){
if (move(b1)){
pro_para=p_para(b0-1,_);
for (b2=b1-1;b2>=0;--b2){
pro_para(b2)=p_para(b0,b2);	
}
pro_para(b1)+=rnorm(0.0,sigma(b1));
newloglik(0)=prior_loglik(pro_para);
if (newloglik(0)> -9999999){
loglikallpro=loglik(data11,SI,pro_para,n_inf,n_sus,with_rm,sep1,sep2);
NumericMatrix tempoutput=loglikallpro(0);
NumericMatrix tempoutput2=loglikallpro(1);
loglik1pro=clone(tempoutput);
loglik2pro=clone(tempoutput2);
newloglik(1)=sum(loglik1pro);
newloglik(2)=sum(loglik2pro);
newloglik(0)+=newloglik(1)+newloglik(2);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
else{
accept_pro=0;	
}
if(gen_binom(accept_pro)){
p_para(b0,b1)=pro_para(b1);
loglik1=clone(loglik1pro);
loglik2=clone(loglik2pro);		
temploglik(2)=newloglik(2);
temploglik(1)=newloglik(1);
temploglik(0)=newloglik(0);
acceptrate(b1)*=(b0-1);
acceptrate(b1)+=1;
acceptrate(b1)/=b0;
}
else{
p_para(b0,b1)=p_para(b0-1,b1);
acceptrate(b1)*=(b0-1);
acceptrate(b1)/=b0;
}
}
else {
p_para(b0,b1)=p_para(b0-1,b1);
}
}

LL1(b0,0)=temploglik(0)-prior_loglik(p_para(b0,_));
LL1(b0,1)=temploglik(1);
LL1(b0,2)=temploglik(2);

// move the matirx to another matrix to store the parameter and compute the correlation matrix
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(b0,moveindex)=p_para(b0,b1);
--moveindex;
}	
}

// update the infection time for each member

for (b1=max_member-1;b1>=0;--b1){

// update the individual parameter without changing infection status
List allupdate=all_update(data11,data1,SI,p_para(b0,_),b1,loglik1,loglik2,n_inf,n_sus,with_rm,sep1,sep2);
NumericMatrix tempoutput1=allupdate(0);
data11=clone(tempoutput1);
NumericMatrix tempoutput3=allupdate(1);
loglik1=clone(tempoutput3);
NumericMatrix tempoutput4=allupdate(2);
loglik2=clone(tempoutput4);
updateacceptrate(b0,b1)=allupdate(3);
}

// update the impute likelihood
LL2(b0,1)=sum(loglik1);
LL2(b0,2)=sum(loglik2);
LL2(b0,0)=LL2(b0,1)+LL2(b0,2);

temploglik(0)=LL2(b0,0)+prior_loglik(p_para(b0,_));
temploglik(1)=LL2(b0,1);
temploglik(2)=LL2(b0,2);

// need update in future, since it only records the index cases
if (with_rm==1){
if (b0>=burnin){
if (b0%(thinning)==thinning-1){	
for (b1=data11.nrow()-1;b1>=0;--b1){
rmrecord(recindex,b1)=data11(b1,7);
}
++recindex;
}
}
}

if (b0%1000==0){
Rcout << "Iteration: " << b0 << std::endl;
}

}


return List::create(_[""]=p_para,
_[""]=LL1,
_[""]=rmrecord,
_[""]=acceptrate,
_[""]=updateacceptrate,
_[""]=data11);
} 

