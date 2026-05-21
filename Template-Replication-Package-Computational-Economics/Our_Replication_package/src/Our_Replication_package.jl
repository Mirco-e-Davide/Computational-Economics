module Our_Replication_package
using DifferentialEquations
using Plots
using ADTypes

greet() = print("Hello World! Earthlings!")

#= MAIN PARAMETERS

eta is the wealth held by experts (net worth, different from psi)

a = Productivity of the capital managed by experts
It indicates the output generated per unit of capital when held by experts

a_ = Productivity of the capital managed by households. It is lower than or equal to a, reflecting the assumption that households are less productive

delta = Depreciation rate of the capital when managed by experts

delta_ = Depreciation rate of the capital when managed by households

rho = Discount rate of the experts. It is strictly greater than the households' discount rate, which drives experts to consume and prevents them from accumulating infinite wealth

r = Discount rate of the households. Since households supply perfectly elastic funding, this acts as the risk-free interest rate of the economy

sigma = Exogenous fundamental risk parameter. It measures the magnitude of the aggregate Brownian shocks that hit the future productivity of capital

=#


# this is the function that will return our output
# we set the default parameters, so that if the user wants to modify them, he can
# otherwise he is not forced to do it
function run_all(; a = 0.11, a_ = 0.05, rho = 0.06, r = 0.05, sigma = 0.025, delta = 0.03, delta_ = 0.08)
    
    # this nested function simply calculates phi(actual rate of new capital creation)
    # and iota (optimal investment rate per unit of capital)
    function investment(q)  
        theta= 10.0 #(investment adjustment cost)

        Phi= (q-1.0)/theta
        iota= Phi+theta*(Phi^2)/2.0

        return Phi, iota
    end

    # we initialize the vairables otherwise they wouldn't be visible outside the for loop
    QL = 0.0
    QR = 10.0
    qmax = 0.0  
    # we look for qmax with a bisection method
    for iter in 1:50
        qmax = (QL + QR) / 2.0
        Phi, iota = investment(qmax)
        
        # expected value for experts
        value = (a - iota) / (r + delta - Phi)
        
        if iota > a  # we would need to invest more than what the assets give back
            QR = qmax           # price too large, so we restrict the interval for qmax
        elseif value > qmax
            QL = qmax           # price too low
        elseif r + delta < Phi   # this means the capital grows faster than it depreciates, so this gives infinite value (economically, not mathematically)
            QR = qmax           # that means the price is too high (high price means more investment and therefore more growing)
        else
            QR = qmax
        end
    end

    # we look for q_ (q at eta=0) with a bisection method
    QL = 0.0
    QR = 10.0
    q_ = 0.0  
    
    for iter in 1:50
        q_ = (QL + QR) / 2.0
        Phi, iota = investment(q_)
        
        # at eta=0 all the wealth is in the hands of the households so we use a_ and delta_
        value = (a_ - iota) / (r + delta_ - Phi)
        
        if iota > a_
            QR = q_             # price too large
        elseif value < q_
            QR = q_             # price too large
        else
            QL = q_             
        end
    end
    p = (qmaxparam=qmax, a = 0.11, a_ = 0.05, rho = 0.06, r = 0.05, sigma = 0.025, delta = 0.03, delta_ = 0.08)

    function fnct!(du, u, p, eta) # p is needed because the function we pass to the ODE must have foru arguments
        # u corresponds to f
        # u[1] = theta, u[2] = theta', u[3] = q, u[4] = q'
        # du corresponds to fp, it contains the 
        
        (; a, a_, rho, r, sigma, delta, delta_) = p # to unpack parameters (p.a mapped to a and so on) (not sure if needed)
        
        Phi, iota = investment(u[3])
        
        # we look for psi (fraction of the economy's total physical capital that is allocated to and managed by experts)
        # with a bisection method
        psi_L = eta
        psi_R = min(u[3] / u[4] + eta, 1.0) #the first term is a mathematical limit in order for the amplification not to be 0
        
        psi = 0.0
        sigma_eta_eta = 0.0
        sigma_q = 0.0
        sigma_theta = 0.0
        risk_premium=0.0
        
        for n in 1:50
            psi = (psi_L + psi_R) / 2.0
            
            # Dynamic amplification due to leverage
            amplification = 1.0 - (u[4] / u[3]) * (psi - eta)
            
            # computation following Proposition II.4
            sigma_eta_eta = sigma* (psi - eta) / amplification #if the amplification is 0 this explodes
            sigma_q = sigma_eta_eta* u[4] / u[3]
            sigma_theta = sigma_eta_eta* u[2] / u[1]
            risk_premium = -sigma_theta * (sigma + sigma_q)
            
            household_premium = (a_ - a) / u[3] + delta - delta_ + risk_premium
            #the household premium is just the risk premium of the experts minus something due to their inefficience in using the capital
            if household_premium > 0.0  # households are willing to buy assets, that means psi must go down
                psi_R = psi
            else
                psi_L = psi
            end
        end
        
        
        mu_q = r - (a - iota) / u[3] - Phi + delta - sigma * sigma_q + risk_premium
        
        mu_eta_eta = -(psi - eta) * (sigma + sigma_q) * (sigma + sigma_q + sigma_theta) + 
                     eta * (a - iota) / u[3] + eta * (1.0 - psi) * (delta_ - delta)
        
        # computing the second derivatives of q and theta wrt eta
        qpp = 2.0 * (mu_q * u[3] - u[4] * mu_eta_eta) / (sigma_eta_eta^2)
        thetapp = 2.0 * ((rho - r) * u[1] - u[2] * mu_eta_eta) / (sigma_eta_eta^2)
        
        # update du
        du[1] = u[2]        
        du[2] = thetapp     
        du[3] = u[4]        
        du[4] = qpp         

    end

    function condition!(out, u, eta, integrator)
        current_qmax= integrator.p.qmaxparam # i call it qmaxparam to distinguish it from qmax found before
        # in out we put the variables to check
        out[1] = current_qmax - u[3]  # if q is higher than qmax
        out[2] = u[2]         # if theta'=0
        out[3] = u[4]         # if q'=0
    end

    # affect! says what to do if a condition above is verified
    last_event= Ref(0) # to capture what condition triggered
    function affect!(integrator, idx)
        if idx isa AbstractArray  # if more conditions verify simultaneously julia returns an array
            e_idx = findfirst(!=(0), idx) # we just take the first
            last_event[] = isnothing(e_idx) ? 0 : e_idx
        else
            last_event[]= idx
        end
        # we terminate immediately no matter what condition is verified
        println("Integration terminated by condition number: ", idx)
        terminate!(integrator)
    end

    # it continously verify condition (the 3 stands for the number of conditions) and eventually runs affect
    cb = VectorContinuousCallback(condition!, affect!, 3,abstol=1e-10)

    # now we try to guess q'(0) with the 'shooting' method
    QL_shoot = 0.0
    QR_shoot = 1e15
    
    # we initialize sol out of the loop in order to use it afterwards
    #local sol 
    # u0 = [theta(0), theta'(0), q(0), q'(0)] so we are evaluating at eta=0
    u0_start = [1.0, -1e10, q_, 0.0]
    prob = ODEProblem(fnct!, u0_start, (0.0, 1.0), p) # p is needed because it gives the value to the p in fnct!
    for iter in 1:50

        last_event[]=0

        guess_qp = (QL_shoot + QR_shoot) / 2.0

        println("Shooting iterazione: $iter/50 | Provando q'(0) = $guess_qp")
        u0_new=copy(u0_start)
        u0_new[4]=guess_qp
        
        # we define the ode with (0,1) as interval of integration
        #prob = ODEProblem(fnct!, u0, (1e-5, 1.0), p) # p is needed because it gives the value to the p in fnct!

        prob_new = remake(prob, u0 = u0_new) #try to rebuild the problem at every iteration

        # Tsit5() it's julia's algorithm for non-stiff problem (better than ode45).
        # cb is the stop condition from above.
        #sol = solve(prob, Tsit5(), callback=cb, reltol=1e-8, abstol=1e-10, maxiters=1e8, dt=1e-6)
        sol= solve(prob_new, Tsit5(), callback=cb,reltol=1e-8, abstol=1e-10)
        # sol[4, end] is the final value of q'. sol[2, end] is the final value of theta'.

        if last_event[]==3
            # if q'(0) reaches 0 before theta'(0) the initial guess must have been to low
            QL_shoot = guess_qp
        else
            # if theta'(0) reaches 0 before q'(0) the initial guess must have been to high
            QR_shoot = guess_qp
        end
    end
    best_qp = (QL_shoot + QR_shoot) / 2.0
    u0_final = copy(u0_start)
    u0_final[4] = best_qp
    prob_final = remake(prob, u0 = u0_final)
    
    # We save this to a new local variable 'sol' which lives in the main function scope
    sol = solve(prob_final, Tsit5(), callback=cb, reltol=1e-8, abstol=1e-10)
    etasol= sol.t
    thetasol= sol[1,:]
    thetapsol= sol[2,:]
    qsol= sol[3,:]
    qpsol=sol[4,:]
    println("eta   = ", etasol[end])
    println("theta = ", thetasol[end])
    println("theta'= ", thetapsol[end])
    println("q     = ", qsol[end])
    println("q'    = ", qpsol[end])
    println("Reason for termination (idx): ", last_event[])

    # this function does the same job as funct! but it also returns dyn, which we need for the plotting.
    # the reason why we couldn't return dyn in fnct! (as in the original code) is that the functions use by julia to solve the ODE don't accept two output 
    # from funct!. Moreover, it would be very heavy in computational cost to calculate dyn at every iteration of the 'shooting' method
    function calc_dynamics(u, p, eta)
        (; a, a_, rho, r, sigma, delta, delta_) = p # to unpack parameters (p.a mapped to a and so on) (not sure if needed)
        theta   = u[1]
        theta_p = u[2]
        q       = u[3]
        q_p     = u[4]
        
        Phi, iota = investment(q)
        
        psi_L = eta
        psi_R = min(q / q_p + eta, 1.0)
        
        psi = 0.0
        amplification = 0.0
        sigma_eta_eta = 0.0
        sigma_q = 0.0
        sigma_theta = 0.0
        risk_premium = 0.0
        household_premium = 0.0
        
        for n in 1:50
            psi = (psi_L + psi_R) / 2.0
            
            amplification = 1.0 - (q_p / q) * (psi - eta)
            
            sigma_eta_eta = sigma * (psi - eta) / amplification
            sigma_q = sigma_eta_eta * (q_p / q)
            sigma_theta = sigma_eta_eta * (theta_p / theta)
            risk_premium = -sigma_theta * (sigma + sigma_q)
            
            household_premium = (a_ - a) / q + delta - delta_ + risk_premium
            
            if household_premium > 0.0 
                psi_R = psi
            else
                psi_L = psi
            end
        end
        
        mu_q = r - (a - iota) / q - Phi + delta - sigma * sigma_q + risk_premium
        
        mu_eta_eta = -(psi - eta) * (sigma + sigma_q) * (sigma + sigma_q + sigma_theta) + 
                     eta * (a - iota) / q + 
                     eta * (1.0 - psi) * (delta_ - delta)
        
        # these lines weren't present in fnct!
        leverage = psi / eta  #(experts' leverage ratio)
        rk = r + risk_premium #(total expected return earned by experts from holding capital)
        r_k = r + household_premium #(total expected return that the less productive households earn if they hold capital)

        dyn = [psi, sigma_eta_eta, sigma_q, mu_eta_eta, mu_q, iota, leverage, rk, r_k]
        
        return dyn
    end

    # PLOTTING
    etaout = sol.t
    N = length(etaout)
    
    # In MATLAB: normalization = fout(N,1); fout(:,1:2) = fout(:,1:2)/normalization;
    # we want theta(eta*)=1
    normalization = sol[1, end] 
    theta_out   = thetasol ./ normalization
    theta_p_out = thetapsol ./ normalization
    q_out       = qsol
    q_p_out     = qpsol
    
    dynout = zeros(N, 9)
    
    for n in 1:N
        u_val = [theta_out[n], theta_p_out[n], q_out[n], q_p_out[n]]
        dynout[n, :] = calc_dynamics(u_val, p, etaout[n])
    end
    
    c = :red # like in original paper
    
    p1 = plot(etaout, q_out, xlabel="\\eta", ylabel="q", color=c, legend=false)
    p2 = plot(etaout, theta_out, xlabel="\\eta", ylabel="\\theta", ylims=(0, 10), color=c, legend=false)
    p3 = plot(etaout[1:N-1], dynout[1:N-1, 1], xlabel="\\eta", ylabel="\\psi", color=c, legend=false)
    p4 = plot(etaout[1:N-1], dynout[1:N-1, 4], xlabel="\\eta", ylabel="\\eta \\mu^\\eta", color=c, legend=false)
    p5 = plot(etaout[1:N-1], dynout[1:N-1, 2], xlabel="\\eta", ylabel="\\eta \\sigma^\\eta", color=c, legend=false)
    p6 = plot(etaout[1:N-1], dynout[1:N-1, 3], xlabel="\\eta", ylabel="\\sigma^q", color=c, legend=false)
    p7 = plot(etaout[1:N-1], dynout[1:N-1, 6], xlabel="\\eta", ylabel="\\iota", color=c, legend=false)
    p8 = plot(etaout[1:N-1], dynout[1:N-1, 7], xlabel="\\eta", ylabel="expert leverage", ylims=(0, 10), color=c, legend=false)
    p9 = plot(etaout[1:N-1], dynout[1:N-1, 8], xlabel="\\eta", ylabel="returns", color=c, label="r^k")
    plot!(p9, etaout[1:N-1], dynout[1:N-1, 9], color=c, label="r_k")
    hline!(p9, [r], color=c, linestyle=:dot, label="r")
    
    # unite the 9 subplots in one single figure
    fig1 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3, 3), size=(900, 700), title="Figure 1")
    
    EXPERT = etaout .* q_out .* theta_out
    HOUSEHOLD = (1.0 .- etaout) .* q_out
    
    fig2 = plot(EXPERT, HOUSEHOLD, xlabel="expert utility", ylabel="household utility", color=c, legend=false, title="Figure 2")
    # this is the output
    return (figure1 = fig1, figure2 = fig2, dynamics = dynout, states = sol)

end

export greet, run_all

end # module Our_Replication_package
