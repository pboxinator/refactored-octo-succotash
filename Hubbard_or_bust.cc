#include<iostream>
using namespace std;
#include "itensor/all.h"
using namespace itensor;


float V(float k){
	// Coulomb interaction.
	if (k!=0) return 1/abs(k);
	else return 0;}


/*-----MAIN-----*/
int main(){
	
bool infinite = false;
int Lam = 5; // Truncation length of interaction.
int Nuc = 1; // Number of unit cells.
int Nx  = 2*Nuc; // iDMRG uses 2 unit cells in the central region.
auto sites = Hubbard(Nx);

auto Hupup = IQMPO(sites);
auto Hdndn = IQMPO(sites);

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
	// Hupup:
	auto& Wupup = Hupup.Aref(n); // Access H-elements by reference using temp W.
	Wupup = IQTensor(dag(sites(n)),prime(sites(n)),row,col);
	Wupup += sites.op("Id",n)  * row(1)     * col(1); // Start node of finite state machine.
	Wupup += sites.op("Id",n)  * row(2+Lam) * col(2+Lam);
	Wupup += sites.op("Nup",n) * row(2)     * col(1);
	for (int m=3; m<2+Lam; ++m){
		Wupup += sites.op("Id",n)  * row(m)     * col(m-1);}
	for (int m=2; m<2+Lam; ++m){
		Wupup += sites.op("Nup",n) * row(2+Lam) * col(m) * V(m-1);}		
	// Hdndn:
	auto& Wdndn = Hdndn.Aref(n);
	Wdndn = IQTensor(dag(sites(n)),prime(sites(n)),row,col);
	Wdndn += sites.op("Id",n)  * row(1)     * col(1);
	Wdndn += sites.op("Id",n)  * row(2+Lam) * col(2+Lam);
	Wdndn += sites.op("Ndn",n) * row(2)     * col(1);
	for (int m=3; m<2+Lam; ++m){
		Wdndn += sites.op("Id",n)  * row(m)     * col(m-1);}
	for (int m=2; m<2+Lam; ++m){
		Wdndn += sites.op("Ndn",n) * row(2+Lam) * col(m) * V(m-1);}		
}


Print(Hupup);
Print(Hdndn);
println("---------END Hupup, Hdndn----------");

//std::vector<IQMPO> Hvec {Hupup,Hdndn};
//auto H = sum(Hvec);
auto H = sum(Hupup,Hdndn);
	
return 0;}
