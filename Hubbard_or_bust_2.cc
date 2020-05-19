#include<iostream>
using namespace std;
#include "itensor/all.h"
using namespace itensor;
#include<iomanip>
#include "directsum.cc"

const bool use_add = true; // Use ITensor's add() function. Else use custom pAdd() --> see below.

const bool infinite = true; // Doesn't work if use_add==true. Does add() not work with iDMRG?
const int Lam = 31; // Truncation length of interaction.
const int Nuc = 121; // Number of unit cells. // iDMRG requires Nx=2*Nuc.
const int Nx = 2*Nuc;

void pAdd(IQMPO& Ha, IQMPO& Hb, IQMPO& Hab, std::vector<IQIndex>& links_a, std::vector<IQIndex>& links_b, std::vector<IQIndex>& links_ab){
/*
	Adds two MPO's Ha + Hb and makes Hab equal to their sum.
	Likewise, the QNLinks links_a and links_b are combined into links_ab.
	Modifies Hab and links_ab by reference.
	Requires "directsum.cc" by Matt Fishman:
		https://gist.github.com/mtfishman/1804a3a473112cb7cb80722eeac47a8b
 
	*Hab is required to be of the same siteSet: IQMPO H = IQMPO(sites);
*/
	std::vector<IQIndex>  newLinks(Nx+1);
	std::vector<IQTensor> newW(Nx+2); // Nx+2 because infinite case stores LH and RH at 0 and Nx+1.
	IQIndexSet linksTemp; // Holds two elements for storing output of directSum().

if (infinite){
	newW.at(1) = directSum(Ha.A(1),Hb.A(1),{links_a.at(0),links_a.at(1)},{links_b.at(0),links_b.at(1)},linksTemp,{"IndexName=","links"});
	newLinks.at(0) = linksTemp[0].dag(); // .dag() makes index direction OUT.
	newLinks.at(1) = linksTemp[1];       // This should be OUT by default.
	Hab.Aref(1) = newW.at(1); // Is this an EXPENSIVE action?
for (int n=2;n<=Nx-1;++n){
	// linksTemp rewritten by reference in the next line:
	newW.at(n) = directSum(Ha.A(n),Hb.A(n),{links_a.at(n-1),links_a.at(n)},{links_b.at(n-1),links_b.at(n)},linksTemp,{"IndexName=","links"});
	newLinks.at(n) = linksTemp[1]; // .dag() makes index direction OUT.
	newW.at(n) = newW.at(n)*delta(linksTemp[0].dag(),newLinks.at(n-1).dag());
	Hab.Aref(n) = newW.at(n);}
	// For n=Nx:
	newW.at(Nx) = directSum(Ha.A(Nx),Hb.A(Nx),{links_a.at(Nx-1),links_a.at(0)},{links_b.at(Nx-1),links_b.at(0)},linksTemp,{"IndexName=","links"});
	newW.at(Nx) = newW.at(Nx)*delta(linksTemp[0].dag(),newLinks.at(Nx-1).dag())*delta(linksTemp[1].dag(),newLinks.at(0));
	Hab.Aref(Nx) = newW.at(Nx);
	// LH (aka Hf):
	newW.at(0) = directSum(Ha.A(0),Hb.A(0),{links_a.at(0)},{links_b.at(0)},linksTemp,{"IndexName=","links"});
	newW.at(0) = newW.at(0)*delta(linksTemp[0].dag(),newLinks.at(0)); // LH has one bond at 0.
	Hab.Aref(0) = newW.at(0);
	// RH (aka Hi)
	newW.at(Nx+1) = directSum(Ha.A(Nx+1),Hb.A(Nx+1),{links_a.at(0)},{links_b.at(0)},linksTemp,{"IndexName=","links"});
	newW.at(Nx+1) = newW.at(Nx+1)*delta(linksTemp[0].dag(),newLinks.at(0).dag()); // RH has one bond at 0 (for intinite).
	Hab.Aref(Nx+1) = newW.at(Nx+1);
}//infinite

else /*finite*/ {
	newW.at(1) = directSum(Ha.A(1),Hb.A(1),{links_a.at(1)},{links_b.at(1)},linksTemp,{"IndexName=","links"});
	newLinks.at(1) = linksTemp[0]; // This should be OUT by default.
	Hab.Aref(1) = newW.at(1);
for (int n=2;n<=Nx-1;++n){
	// linksTemp rewritten by reference in the next line:
	newW.at(n) = directSum(Ha.A(n),Hb.A(n),{links_a.at(n-1),links_a.at(n)},{links_b.at(n-1),links_b.at(n)},linksTemp,{"IndexName=","links"});
	newLinks.at(n) = linksTemp[1]; // .dag() makes index direction OUT.
	newW.at(n) = newW.at(n)*delta(linksTemp[0].dag(),newLinks.at(n-1).dag());
	Hab.Aref(n) = newW.at(n);}
	// For n=Nx:
	newW.at(Nx) = directSum(Ha.A(Nx),Hb.A(Nx),{links_a.at(Nx-1)},{links_b.at(Nx-1)},linksTemp,{"IndexName=","links"});
	newW.at(Nx) = newW.at(Nx)*delta(linksTemp[0].dag(),newLinks.at(Nx-1).dag());
	Hab.Aref(Nx) = newW.at(Nx);
}//finite

	// Written by reference:
	links_ab = newLinks;
return;}

float V(float k){
	// Coulomb interaction.
	if (k!=0) return 1/abs(k);
	else return 0;}


int main(){
//////////////////////////////////////////////
/////////// Define links and make MPO ////////
auto sites = Hubbard(Nx);
auto Ha = IQMPO(sites);
auto Hb = IQMPO(sites);
// link{q0}:
std::vector<IQIndex> links(Nx+1);
std::vector<Index> q0(Nx+1);
for (int l=0; l<=Nx; ++l){	// Links carry difference of physical site QN's.
	q0.at(l) = Index(nameint("q0_",l),2+Lam);
	links.at(l) = IQIndex(nameint("l",1),
							q0[l],QN("Sz=", 0,"Nf=",0),
							Out); // link{q0}
}
// If infinite: Store final link at 0. If finite: make the final link Nx.
IQIndex const& last = (infinite ? links.at(0) : links.at(Nx));
for (int n=1; n<=Nx; ++n){
// link{q0}:
	auto row = dag(links.at(n-1));
	auto col = (n==Nx ? last : links.at(n));
	// Ha:
	auto& Wa = Ha.Aref(n); // Access H-elements by reference using temp W.
	Wa = IQTensor(dag(sites(n)),prime(sites(n)),row,col);
	Wa += sites.op("Id",n)  * row(1)     * col(1); // Start node of finite state machine.
	Wa += sites.op("Id",n)  * row(2+Lam) * col(2+Lam);
	Wa += sites.op("Nup",n) * row(2)     * col(1);
	for (int m=3; m<2+Lam; ++m){
		Wa += sites.op("Id",n)  * row(m)     * col(m-1);}
	for (int m=2; m<2+Lam; ++m){
		Wa += sites.op("Nup",n) * row(2+Lam) * col(m) * V(m-1);}		
	// Hb:
	auto& Wb = Hb.Aref(n);
	Wb = IQTensor(dag(sites(n)),prime(sites(n)),row,col);
	Wb += sites.op("Id",n)  * row(1)     * col(1);
	Wb += sites.op("Id",n)  * row(2+Lam) * col(2+Lam);
	Wb += sites.op("Ndn",n) * row(2)     * col(1);
	for (int m=3; m<2+Lam; ++m){
		Wb += sites.op("Id",n)  * row(m)     * col(m-1);}
	for (int m=2; m<2+Lam; ++m){
		Wb += sites.op("Ndn",n) * row(2+Lam) * col(m) * V(m-1);}		
}
// Edge vectors:
auto LH = setElt(links.at(0)(2+Lam));
auto RH = setElt(dag(last)(1));
if (not infinite)
{
	// Multiply first and last MPO tensors by respective edge vectors:
	Ha.Aref(1)  *= LH;
	Ha.Aref(Nx) *= RH;
	Hb.Aref(1)  *= LH;
	Hb.Aref(Nx) *= RH;
}
else
{
	// Store edge vectors just before and after last sites:
	Ha.Aref(0)    = LH;
	Ha.Aref(Nx+1) = RH;
	Hb.Aref(0)    = LH;
	Hb.Aref(Nx+1) = RH;
}


//Print(Ha);
//Print(Hb);
//println("---------END Ha & Hb----------");


//////////////////////////////////////////////
/////////// add MPOs /////////////////////////
IQMPO H = IQMPO(sites);
if (use_add){cout << "Using ITensor::add...\n"; 
std::vector<IQMPO> Hvec {Ha,Ha};
H = sum(Hvec);} // Note: Don't define as "IQMPO H = sum(Hvec)" because this gives a bug.
else {cout << "Using pAdd...\n";
pAdd(Ha,Hb,H,links,links,links);}



//////////////////////////////////////////////
/////////// DMRG /////////////////////////////
auto sweeps = Sweeps(21);
sweeps.maxm() = 20,80,140,200;
sweeps.cutoff() = 1E-10,Args("Repeat",10),1E-14;
sweeps.niter() = 3,2,5,5;

auto state = InitState(sites);
for (int i=1; i<=Nx; ++i)
{
	if (i%2 == 1)
		state.set(i,"UpDn");
	else
		state.set(i,"Emp");
}
auto psi = IQMPS(state);

if (!infinite){
auto Hp_ampo = AutoMPO(sites);
for (int n=1; n<Nx; ++n){
for (int k=1; k<=Lam; ++k){
	if (n+k > Nx) continue;
	Hp_ampo += V(k),"Nup",n,"Nup",n+k;
	Hp_ampo += V(k),"Ndn",n,"Ndn",n+k;
	}} auto Hp = IQMPO(Hp_ampo);

auto en  = dmrg(psi,H ,sweeps); // Uses add() or pAdd() to combine Ha and Hb into H.
auto enp = dmrg(psi,Hp,sweeps); // Uses autoMPO to construct H.
cout << "energy per site (autoMPO) = "        << std::fixed << std::setprecision(12) << enp/Nx << endl;
cout << "energy per site (using add/pAdd) = " << std::fixed << std::setprecision(12) << en/Nx  << endl;}

else {
auto res = idmrg(psi,H,sweeps,{"OutputLevel",1});
printfln("\nGround state energy / site = %.20f",res.energy/Nx);}
	
return 0;}
