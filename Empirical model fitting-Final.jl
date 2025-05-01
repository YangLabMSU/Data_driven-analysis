using RHEOS
using PyPlot
using Statistics
using CurveFit
using CSV
using DataFrames, XLSX


cd(raw"#The location of the input files#")

#The time point file inclues the point afterwhich the strain is kept constant (t_hold) and the transition point (t_trans)
#after which the stress will start to increase.

T =  CSV.read("#File containing the t_hold and t_trans for all the curves#.csv", DataFrame);
n = convert(Int64,size(T,2)/2); #The number of curves
ϵ0 = 0.5 #The applied strain

result1 = zeros(5,n);
result2 = zeros(5,n);
result3 = zeros(5,n);

# n is the number of experiments.

for i = 1:n

    data = RHEOS.importcsv("#File containing time, strain and stress#.csv", t_col = 3*i-2, ϵ_col = 3*i-1, σ_col = 3*i);
    data = resample(data);


    fig, ax = subplots(2, 3, figsize = (30, 20))
    ax[1, 1].plot(data.t, data.σ, ".", color = "green")
    ax[1, 1].set_title("Stress-Time", fontsize = 40)
    ax[1, 2].plot(data.t, data.ϵ)
    ax[1, 2].set_title("Strain-time", fontsize = 40)

    ## Separating the curve into two pieces (Relaxed and Tensioned). The relaxed part
    # is used to calculate initial values for A, B, α, and τ. The tensioned part
    # is used to calculate initial value for C.

    #Relaxed part
    t_h = T[1,2*i-1]; # The hold point

    T2 = findmin(data.σ[t_h:end])


    t_trans = T2[2]+t_h-1; # The transition point


    time_on1 = data.t[t_h]; # First time point after hold
    time_off1 = data.t[t_trans]; # Time at the transition point
    rel = cutting(data::RheoTimeData, time_on1::Real, time_off1::Real);

    #x_rel is a vector including the time points for relaxed part of the curve starting at time zero.

    x_rel = rel.t.- data.t[t_h];

    # Initial value for B

    b = minimum(data.σ[t_h:end]);



    t_b = findall(x-> x==b, data.σ); #The point where the stress has its minimum value
    T_b = data.t[t_trans]; #Time at the initial value for B

    # Initial value for α

    t0 = t_h;
    t1 = t_trans-2;
    fit1 = CurveFit.curve_fit(LinearFit, log.(data.t[t0:t1]), log.(data.σ[t0:t1].- b));
    α_i = fit1.(1) - fit1(0);
    y1 = fit1.(log.(data.t[t0:t1]));
    ax[1, 3].plot(log.(data.t[t0:t1]),log.(data.σ[t0:t1] .- b), "o", log.(data.t[t0:t1]), y1,"r-", linewidth=3)
    ax[1, 3].set_title("α", fontsize = 40)

    if -α_i>1
        α_i= -0.9;
    end


    #Cutting the curve after the point σ(min)+σ(min)*0.4

    t3 = maximum(findall(data.σ .< 1.4.*data.σ[t_b]))


    time_off3 = data.t[t3];

    time_on3 = data.t[1];

    Data_Cut = cutting(data::RheoTimeData, time_on3::Real, time_off3::Real);

    # Initial value for A. (I used the initial values calculated for α, τ and b. Then, fitted the relaxed part of
    #the curve to calculate the initial value for A.)

    aa = data.σ[t_h] - b

    using LsqFit
    @. model(t, p) = p[1]*t.^(α_i)+b;
    P0 = [aa]
    fit = LsqFit.curve_fit(model, x_rel, rel.σ, P0)
    a = fit.param
    y0 = a.*(x_rel).^(α_i).+b;


    # Tensioned part of the curve is deined between the transition point and the point where σ(min)+σ(min)*0.4

    time_on2 = data.t[t_trans-1];

    TEN = cutting(data::RheoTimeData, time_on2::Real, time_off3::Real);

    # Initial value for C

    fit3 = CurveFit.curve_fit(LinearFit, TEN.t, (TEN.σ.- b));
    c = fit3.(1) - fit3(0);
    y3 = fit3.(TEN.t)

    ax[2,1].plot(TEN.t, TEN.σ.- b, "o", TEN.t, y3, "r-", linewidth = 3)
    ax[2,1].set_title("C", fontsize = 40)


    if c<=0
        c = 0;
    end

    # Reporting the initial values for the parameters

    result1[1,i] = a[1]/ϵ0; result1[2,i]=aa/ϵ0; result1[3,i]=-α_i; result1[4,i] = b/ϵ0;result1[5,i] =c/ϵ0 ;

    ## Fitting the curve with the convolution integral

    # Lower bounds
    Lo = (A = 0, α = 0, B = 0.85 * b/ϵ0 , C = 0)# B was allowed to vary by 15% to optimize the fits.
    # Upper bounds
    Hi = (A = Inf, α = 1, B = 1.15 * b/ϵ0, C = Inf) # B was allowed to vary by 15% to optimize the fits.
    # Initial parameters
    P0 = (A = a[1]/ϵ0, α = -α_i, B = b/ϵ0, C = c/ϵ0)

    PowerLawEmpirical = RheoModelClass(
    # Model name
    name = "power_empirical",
    # Model parameters,
    p = (:A, :α, :B, :C),
    #Relaxation modulus
    G = quote
        (A * t^(-α) + B + C * t)
    end,
    # Network
    info = "Empirical model",
    )

    # Optimizer LN_BOBYQA

    PLE_model1 =
    modelfit(Data_Cut, PowerLawEmpirical, strain_imposed, lo = Lo, hi = Hi, p0 = P0, optmethod=:LN_BOBYQA)
    PowerLawEmpirical_predict1 = extract(Data_Cut, strain_only)
    PowerLawEmpirical_predict1 = modelpredict(PowerLawEmpirical_predict1, PLE_model1)

    PP1 = PLE_model1.fixedparams;
    result2[1,i] = PP1.A; result2[2,i]=PP1.α; result2[3,i]=PP1.B; result2[4,i] = PP1.C;
    #result2[5,i] = PP.C*((PP.τ)^(1+PP.α))/PP.A;

    Error1 = Data_Cut.log[end].info.error;
    result2[5,i] = Error1;

    n1 = convert(Int64,size(Data_Cut.t,1))

    result3 = zeros(n1,8)

    result3[:,1] = PowerLawEmpirical_predict1.t;
    result3[:,2] = PowerLawEmpirical_predict1.σ;

 


    ax[2,2].plot(Data_Cut.t, Data_Cut.σ, ".", color = "green")
    ax[2,2].plot(PowerLawEmpirical_predict1.t, PowerLawEmpirical_predict1.σ, ".", color = "red");#, markersize = 20);
    ax[2,2].set_title("Experimental loading condition", fontsize = 35)
    ax[2,2].legend(["Experiment", "Empirical model"], fontsize = 20)



    gcf()


    savefig(string(i)*"_100nmps_NEWV2.png")
    plt.clf()

    plot(data.t, data.σ, ".", color = "green", markersize = 20)
    plot(PowerLawEmpirical_predict1.t, PowerLawEmpirical_predict1.σ, ".", color = "red", markersize = 25);
    title("Experimental loading condition", fontsize = 35)
    legend(["Experiment", "Empirical model"], fontsize = 30)
    savefig(string(i)*"_Drug studies_BOB_NEWV2.png")
    plt.clf()

end

RESULT1 = DataFrames.DataFrame(result1, :auto);
RESULT2 = DataFrames.DataFrame(result2, :auto);
RESULT3 = DataFrames.DataFrame(result3, :auto);
XLSX.writetable("Test.xlsx", Positive=(collect(DataFrames.eachcol(RESULT2)), DataFrames.names(RESULT2)),Initial_Values=(collect(DataFrames.eachcol(RESULT1)), DataFrames.names(RESULT1)))#, Fittings=(collect(DataFrames.eachcol(RESULT6)), DataFrames.names(RESULT6)))
XLSX.writetable("Test-fit.xlsx", Fittimg=(collect(DataFrames.eachcol(RESULT3)), DataFrames.names(RESULT3)))