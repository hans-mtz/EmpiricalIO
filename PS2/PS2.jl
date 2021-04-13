## Set up
using DataFrames
using JuMP
using GLM
using Gurobi
using CSV

dir="/Volumes/SSD Hans/Github/EmpiricalIO"

## Reading Data

df = CSV.read(dir*"/PS2/df.csv", DataFrame; missingstring = "NA")

# generating intermediate inputs

df.int = df.rmats + df.renergy + df.rserv

df[:,:int]

## Estimating base ACF model

 # first stage

 
