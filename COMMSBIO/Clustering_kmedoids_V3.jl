using CSV
using Clustering
using Statistics
using PyPlot
import Plots
using DataFrames, XLSX

cd(raw"Location of the input files")


function load_all()
  d=CSV.File("#File containing the value of B in the first column and C in the second column.csv#", header=0)
  d1=filter(!ismissing,d.Column1)
  d2=filter(!ismissing, d.Column2)
  d=hcat(d1,d2)
  return(d)
end


function semilog_norm(d)
  d2 = copy(d)
  B = d2[:,1]
  logC = log.(d2[:,2])
  # Normalise data
  B = B / std(B)
  logC = logC / std(logC)
  d2[:,1] = B
  d2[:,2] = logC

  return(d2)
end

function dmatrix(d)
  l = length(d[:,1])
  m=zeros(l,l)
  for i in 1:l
      for j in 1:l
         m[i,j] = sqrt( (d[j,1]-d[i,1])^2 + (d[j,2]-d[i,2])^2 )
      end
  end
  return(m)
end



function kmedoids_analysis(d,k)
  m = dmatrix(d)
  r=kmedoids(m,k,init=:kmcen)
  return(r.assignments)
end




function plot_clusters(d,c)
  (main_fig, main_ax) = subplots(1, 1, figsize=(12,9))
  scatter(d[:,2],d[:,1],c=r)
  xlim((-10, 2))
  ylim((-6, 9))
end



d=load_all()
dd = d[vec(mapslices(col -> any(col[2] .!= 0.0), d, dims = 2)), :] #Removing the rows that contain C=0
dn = semilog_norm(dd)

n = size(dd,1); #Number of elements in each column
result1 = zeros(n,1);

# kmedoids analysis

r=kmedoids_analysis(dn,2)
plot_clusters(dn,r)
savefig("Test-Kmedoids.png")
plt.clf()

result1 = r;

#Export the results of Clustering
RESULT1 = DataFrames.DataFrame([result1], :auto);
XLSX.writetable("Test-Kmedoids.xlsx", kmeans=(collect(DataFrames.eachcol(RESULT1)), DataFrames.names(RESULT1)));
