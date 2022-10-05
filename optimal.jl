using PATHSolver, XLSX, Complementarity, JuMP
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

discr=0.5
dr=1/discr
adr=9.306414218271248
RPS=0.3
CAP=1.0

## sets
regions=["BJ", "TJ", "HB", "SX", "NM", "LN", "JL", "HLJ", "SH", "JS", "ZJ", "AH", "FJ", "JX", "SD", "HN", "HUB", "HUN", "GD", "GX", "HIN", "CQ", "SC", "GZ", "YN", "XZ", "SHX", "GS", "QH", "NX", "XJ"]

pg_techs =["C","NG","NU","HD","WD","CPV","DPV","BM"]

## parameters
GEF=Dict{String, Float64}("C"=>1.074000,"NG"=>0.354000,"NU"=>0.000000,"HD"=>0.000000,"WD"=>0.000000,"CPV"=>0.000000,"DPV"=>0.000000,"BM"=>0.000000)
CAPCOST=Dict{String, Float64}("C"=>0.254351,"NG"=>0.151787,"NU"=>0.841226,"HD"=>0.559336,"WD"=>0.589898,"CPV"=>0.554667,"DPV"=>0.504841,"BM"=>0.854851)
GVC=Dict{String, Float64}("C"=>0.000156,"NG"=>0.000442,"NU"=>0.000122,"HD"=>0.000000,"WD"=>0.000000,"CPV"=>0.000000,"DPV"=>0.000000,"BM"=>0.000534)
CG=Dict{String, Int}("C"=>0,"NG"=>0,"NU"=>0,"HD"=>0,"WD"=>1,"CPV"=>1,"DPV"=>1,"BM"=>1)
DEM=Dict{String, Float64}("BJ"=>155375.635595,"TJ"=>117164.855413,"HB"=>457849.309853,"SX"=>293854.691058,"NM"=>469503.194858,"LN"=>334492.480613,"JL"=>109043.694427,"HLJ"=>110850.992292,"SH"=>216413.745135,"JS"=>846498.260763,"ZJ"=>626143.491577,"AH"=>294930.326633,"FJ"=>319619.398785,"JX"=>206731.790411,"SD"=>840403.694303,"HN"=>494511.436727,"HUB"=>299719.641798,"HUN"=>252522.512299,"GD"=>875149.312833,"GX"=>235659.973341,"HIN"=>45226.231736,"CQ"=>161187.045163,"SC"=>355797.974916,"GZ"=>205124.862539,"YN"=>232469.861824,"XZ"=>10357.022993,"SHX"=>207806.718655,"GS"=>193503.162705,"QH"=>110794.035883,"NX"=>140327.007416,"XJ"=>320874.137591)

initial_data=XLSX.readxlsx("./ee2020data.xlsx")
powercapacity=initial_data["capacity!A1:I32"]
uh=initial_data["uh!A1:I32"]
tcap=initial_data["tcap!A1:AF32"]
distance=initial_data["distance!A1:AF32"]
totalcap=initial_data["total!A1:I32"]

GIC=Dict{Tuple, Float64}()
GUH=Dict{Tuple, Float64}()

TCAP=Dict{Tuple, Float64}()
TEF=Dict{Tuple, Float64}()
GRIDCOST=Dict{Tuple, Float64}()
TOTAL_PG_CAP=Dict{Tuple, Float64}()

for i in 2:32
	for j in 2:32
		r=tcap[i,1]
		r1=tcap[1,j]
		cap=tcap[i,j]
		dist=distance[i,j]+distance[j,i]
		if r==r1
			cap=1000
			dist=0
		else
			cap=max(tcap[i,j],tcap[j,i],tcap[i,j]+tcap[j,i])
		end
		effi=1-dist/1000*0.03
		gcost=discr/(1+discr)/(1-(1+discr)^(-40))*dist*0.786/1000*1.03/5000
		TCAP[r,r1]=cap
		TEF[r,r1]=effi
		GRIDCOST[r,r1]=gcost
	end

	for j in 2:9
		r=powercapacity[i,1]
		t=powercapacity[1,j]
		v=powercapacity[i,j]
		if ismissing(v)
			v=0
		end
		GIC[r,t]=v

		u=uh[i,j]
		if ismissing(u)
			u=0
		end
		GUH[r,t]=u

		tot=totalcap[i,j]
		TOTAL_PG_CAP[r,t]=tot
	end
end

ee2020_model=MCPModel()
@variable(ee2020_model, 0<= newpg[pg in pg_techs, r in regions] <= TOTAL_PG_CAP[r,pg])
@variable(ee2020_model, pg_s[pg in pg_techs,r in regions,r1 in regions] >= 0)
@variable(ee2020_model, pg_b[pg in pg_techs,r in regions,r1 in regions] >= 0)
@variable(ee2020_model, rec_s[r in regions] >= 0)
@variable(ee2020_model, rec_b[r in regions] >= 0)
@variable(ee2020_model, caq_s[r in regions] >= 0)
@variable(ee2020_model, caq_b[r in regions] >= 0)
@variable(ee2020_model, pri_pg[pg in pg_techs,r in regions,r1 in regions])
@variable(ee2020_model, pri_rec)
@variable(ee2020_model, pri_caq)

@variable(ee2020_model, a_cgoutputmax[pg in pg_techs,r in regions] <= 0)
@variable(ee2020_model, a_selllimit[r in regions,r1 in regions] <= 0)
@variable(ee2020_model, a_capmeet[r in regions] <= 0)
@variable(ee2020_model, a_caqselllimit[r in regions] <= 0)

@variable(ee2020_model, b_demandmeet[r in regions])
@variable(ee2020_model, b_rpsmeet[r in regions] <= 0)
@variable(ee2020_model, b_recselllimit[r in regions] <= 0)

@mapping(ee2020_model, cgoutputmax[pg in pg_techs,r in regions], sum(pg_s[pg,r,r1] for r1 in regions)-(GIC[r,pg]+newpg[pg,r])*GUH[r,pg])
@mapping(ee2020_model, selllimit[r in regions,r1 in regions], sum(pg_s[pg,r,r1] for pg in pg_techs)-TCAP[r,r1]*8000)
@mapping(ee2020_model, capmeet[r in regions],caq_s[r]-caq_b[r]-sum(pg_s[pg,r,r1]*(CAP-GEF[pg]) for pg in pg_techs for r1 in regions))
@mapping(ee2020_model, caqselllimit[r in regions],caq_s[r]-sum(pg_s[pg,r,r1]*CAP for pg in pg_techs for r1 in regions))

@mapping(ee2020_model, kkt_newpg[pg in pg_techs,r in regions],dr*CAPCOST[pg]+a_cgoutputmax[pg,r]*GUH[r,pg])
@mapping(ee2020_model, kkt_pg_s[pg in pg_techs,r in regions,r1 in regions],dr*(GVC[pg]-pri_pg[pg,r,r1])-a_cgoutputmax[pg,r]-a_selllimit[r,r1]+a_capmeet[r]*(CAP-GEF[pg])+a_caqselllimit[r]*CAP)
@mapping(ee2020_model, kkt_caq_b[r in regions],dr*pri_caq+a_capmeet[r])
@mapping(ee2020_model, kkt_caq_s[r in regions],-dr*pri_caq-a_capmeet[r]-a_caqselllimit[r])

@mapping(ee2020_model,demandmeet[r in regions], sum(pg_b[pg,r,r1]*TEF[r1,r] for pg in pg_techs for r1 in regions)-DEM[r]*0)
@mapping(ee2020_model,rpsmeet[r in regions], RPS*DEM[r]+rec_s[r]-sum(pg_b[pg,r,r1]*TEF[r1,r]*CG[pg] for pg in pg_techs for r1 in regions)-rec_b[r])
@mapping(ee2020_model,recselllimit[r in regions],rec_s[r]-sum(pg_b[pg,r,r1]*TEF[r1,r]*CG[pg] for pg in pg_techs for r1 in regions))

@mapping(ee2020_model, kkt_pg_b[pg in pg_techs,r in regions,r1 in regions],dr*(pri_pg[pg,r1,r]+GRIDCOST[r1,r])
-b_demandmeet[r]*TEF[r1,r]+b_rpsmeet[r]*TEF[r1,r]*CG[pg]+b_recselllimit[r]*TEF[r1,r]*CG[pg])
@mapping(ee2020_model, kkt_rec_b[r in regions],dr*pri_rec+b_rpsmeet[r])
@mapping(ee2020_model, kkt_rec_s[r in regions],-dr*pri_rec-b_rpsmeet[r]-b_recselllimit[r])

@mapping(ee2020_model, mc_pg[pg in pg_techs,r in regions,r1 in regions],pg_s[pg,r,r1]-pg_b[pg,r1,r])
@mapping(ee2020_model, mc_rec,sum(rec_s[r]-rec_b[r] for r in regions))
@mapping(ee2020_model, mc_caq,sum(caq_s[r]-caq_b[r] for r in regions))

@complementarity(ee2020_model, cgoutputmax, a_cgoutputmax)
@complementarity(ee2020_model, selllimit, a_selllimit)
@complementarity(ee2020_model, capmeet, a_capmeet)
@complementarity(ee2020_model, caqselllimit, a_caqselllimit)

@complementarity(ee2020_model, kkt_newpg,newpg)
@complementarity(ee2020_model, kkt_pg_s,pg_s)
@complementarity(ee2020_model, kkt_caq_b,caq_b)
@complementarity(ee2020_model, kkt_caq_s,caq_s)

@complementarity(ee2020_model, demandmeet,b_demandmeet)
@complementarity(ee2020_model, rpsmeet,b_rpsmeet)
@complementarity(ee2020_model, recselllimit,b_recselllimit)

@complementarity(ee2020_model, kkt_pg_b, pg_b)
@complementarity(ee2020_model, kkt_rec_b, rec_b)
@complementarity(ee2020_model, kkt_rec_b, rec_b)


@complementarity(ee2020_model, mc_pg,pri_pg)
@complementarity(ee2020_model, mc_rec,pri_rec)
@complementarity(ee2020_model, mc_caq,pri_caq)

status = solveMCP(ee2020_model, solver=:PATH, convergence_tolerance=1e-8, output="yes", time_limit=3600)