using RHEOS
using PyPlot
using Statistics
using CurveFit
using CSV
using DataFrames, XLSX


cd(raw"Location of the input files")

#The time point file inclues the point afterwhich the strain is kept constant (t_hold) and the transition point (t_trans) after which the stress will start to increase.

T =  CSV.read("#File containing the value of B in the first column and C in the second column#.csv", DataFrame);
n = convert(Int64,size(T,2)/2); #The number of curves
ϵ0 = 0.5 #The applied strain


result1 = zeros(8,n);
result2 = zeros(8,n);
result3 = zeros(1000,8)


for i = 2 # Data 3 and 4 are 5umps and 10umps, respectively. 
            #Data 1 and 2 are 1umps and 2umps, respectively.

    data = RHEOS.importcsv("#File containing time, strain and stress#.csv", t_col = 3*i-2, ϵ_col = 3*i-1, σ_col = 3*i);
    data = resample(data);


    t_h = T[1,2*i-1]; # The hold point

    T1 = findmin(data.σ[t_h:end])

    t_trans = T1[2]+t_h-1; # The transition point


    # Tensioned part of the curve is defined between the transition point and the end point

    time_on2 = data.t[t_trans-1];

    time_off2 = data.t[end];

    TEN = cutting(data::RheoTimeData, time_on2::Real, time_off2::Real);

    #Cutting the curve after the point σ(min)+σ(min)*0.4
    
    t2 = maximum(findall(data.σ .< 1.4.*data.σ[t_trans]))

    time_off1 = data.t[t2];

    #time_on1 = data.t[1];

    Data_Cut = cutting(data::RheoTimeData, time_on2::Real, time_off1::Real); #The relaxation part is also removed.

    #Initial value for E1 

    e1 = minimum(data.σ[t_h:end]);

 

    ## Fitting the curve with the convolution integral

    # Lower bounds
    Lo = (E1 = 824.6188394, etaA = 546000*4, SA = 957.6142455*2)
    # Upper bounds
    Hi = (E1 = 824.6188394, etaA = 546000*4, SA = 957.6142455*2) 
    # Initial parameters
    P0 = (E1 = 824.6188394, etaA = 546000*4, SA = 957.6142455*2)

    ActiveStress = RheoModelClass(
    # Model name
    name = "active_stress",
    # Model parameters,
    p = (:E1, :etaA, :SA),

    #Relaxation modulus
    G = quote
        l0=0.000020; #Initial length
        A = 0.000000000120; #Cross sectional area
        kb=0.11 # Stiffness of the sensing island
        ϵ0 = 0.5
        KEQ = E1*kb*A/(E1*A+kb*l0)
        TAU = etaA*A/(KEQ*l0)
        ((1/ϵ0)*((((KEQ*l0*ϵ0/A))*(exp(-t/TAU)))))
    end,
    # Network
    info = "Active model",
    )

    #Fitting the model on the tensioned part of the curve
    # Optimizer LN_BOBYQA

    PLE_model1 =
    modelfit(TEN, ActiveStress, strain_imposed, lo = Lo, hi = Hi, p0 = P0, optmethod=:LN_BOBYQA)
    PowerLawEmpirical_predict1 = extract(TEN, strain_only)
    PowerLawEmpirical_predict1 = modelpredict(PowerLawEmpirical_predict1, PLE_model1)

    PP1 = PLE_model1.fixedparams;
    result1[1,i] = PP1.E1; result1[2,i]=PP1.etaA; result1[3,i]=PP1.SA;

    l0=0.000020; #Initial length
    A = 0.000000000120; #Cross sectional area
    kb=0.11 # Stiffness of the sensing island

    keq = PP1.E1*A*kb/(PP1.E1*A+kb*l0); result1[4,i] = keq;
    Tau1 = (PP1.etaA*A)/(keq*l0); result1[5,i] = Tau1/60; #Tau is the characteristic time
    B = keq*l0/A; result1[6,i] = B;
    result1[7,i] = (1/(ϵ0))*(PP1.SA-(keq*l0*ϵ0/A))/Tau1; # C predicted in the active model

    Error1 = TEN.log[end].info.error;
    result1[8,i] = Error1;

    fig, ax = subplots(figsize = (30, 20))

    plot(data.t, data.σ, ".", color = "green", markersize = 20)
    plot(PowerLawEmpirical_predict1.t, PowerLawEmpirical_predict1.σ, ".", color = "red", markersize = 25);
    title("Experimental loading condition", fontsize = 35)
    legend(["Experiment", "Active model fitting"], fontsize = 30)
    savefig(string(i)*"#Name#.png")
    plt.clf()

    #Exporting the fitting curve


    result3[(t_trans-1:end),2*i-1] = PowerLawEmpirical_predict1.t;
    result3[(t_trans-1:end),2*i] = PowerLawEmpirical_predict1.σ;



    #fitting the model on the curve that we cut at the point σ(min)+σ(min)*0.4
    # Optimizer LN_BOBYQA

    PLE_model2 =
    modelfit(Data_Cut, ActiveStress, strain_imposed, lo = Lo, hi = Hi, p0 = P0, optmethod=:LN_BOBYQA)
    PowerLawEmpirical_predict2 = extract(Data_Cut, strain_only)
    PowerLawEmpirical_predict2 = modelpredict(PowerLawEmpirical_predict2, PLE_model2)

    PP2 = PLE_model2.fixedparams;
    result2[1,i] = PP2.E1; result2[2,i]=PP2.etaA; result2[3,i]=PP2.SA;
    keq2 = PP2.E1*A*kb/(PP2.E1*A+kb*l0); result2[4,i] = keq2;
    Tau2 = PP2.etaA*A/(keq2*l0); result2[5,i] = Tau2/60; #Tau is the characteristic time
    B = keq*l0/A; result2[6,i] = B;
    result2[7,i] = (1/(ϵ0))*(PP2.SA-(keq*l0*ϵ0/A))/Tau2; # C predicted in the active model


    Error2 = Data_Cut.log[end].info.error;
    result2[8,i] = Error2;

    fig, ax = subplots(1, 1, figsize = (30, 20))

    plot(data.t, data.σ, ".", color = "green", markersize = 20)
    plot(PowerLawEmpirical_predict2.t, PowerLawEmpirical_predict2.σ, ".", color = "red", markersize = 25);
    title("Experimental loading condition", fontsize = 35)
    legend(["Experiment", "Active model fitting"], fontsize = 30)
    savefig(string(i)*"#Name#.png")
    plt.clf()





end

RESULT1 = DataFrames.DataFrame(result1, :auto);
RESULT2 = DataFrames.DataFrame(result2, :auto);
RESULT3 = DataFrames.DataFrame(result3, :auto);
XLSX.writetable("#Name#.xlsx", TEN=(collect(DataFrames.eachcol(RESULT1)), DataFrames.names(RESULT1)), Cut=(collect(DataFrames.eachcol(RESULT2)), DataFrames.names(RESULT2)))
XLSX.writetable("#Name#.xlsx", Fittimg=(collect(DataFrames.eachcol(RESULT3)), DataFrames.names(RESULT3)))