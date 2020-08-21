using Splines2, GLM, Random
Random.seed!(12345);

x = collect(range(0.0, length=301, stop=2.0*pi));
y = sin.(x)+randn(length(x)); 

ns1 = Splines2.bs_(x,df=5,intercept=false) # this is a function

newx = collect(0.0:0.5:3.5);


d = DataFrames.DataFrame(x=x,y=y);
fit2 = lm(@formula(y~bs2(x,5)+1),d) # equivalent to fit1 with nicer labels
predict(fit2, hcat(ones(8),ns1(newx))) # safe predictions
predict(fit2, DataFrames.DataFrame(x=newx)) # unsafe predictions!
