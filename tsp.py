import numpy as np
import pyomo.opt
from pyomo.environ import *

# Create the (x,y) coordinates for n cities
num_of_cities = 1000
city_coords = np.random.rand(num_of_cities,2)
# D = {}
D = np.empty((num_of_cities,num_of_cities))
for i in range(num_of_cities):
    for j in range(num_of_cities):
        D[i,j] = np.sqrt( (city_coords[i,0]-city_coords[j,0])**2 \
            + (city_coords[i,1]-city_coords[j,1])**2)

print("D is created")
model = ConcreteModel()

model.name = 'TSP'

# model.cities = Set(initialize=RangeSet(num_of_cities), doc = 'Set of cities')
# model.cities = RangeSet(num_of_cities)
model.cities = range(num_of_cities)
# model.d = Param(model.cities, model.cities, initialize=D)

model.x = Var(model.cities,model.cities,
    within=Binary,
    doc='Arc Variables')

# model.slack = Var(model.cities,
#     within=NonNegativeIntegers,
#     doc='Slack Variables')

def in_out_con(model,i):
    total = 0
    for j in model.cities:
        total += model.x[i,j]
    return total == 2

model.in_out_con = ConstraintList()
for i in model.cities:
    model.in_out_con.add(
        sum( model.x[i,j] for j in model.cities ) == 2
    )

def own_city(model,i):
    return model.x[i,i]==0

model.own_city = ConstraintList()
for i in model.cities:
    model.own_city.add(
        model.x[i,i] == 0
    )

def symmetric(model,i,j):
    return model.x[i,j] == model.x[j,i]

model.symmetric = ConstraintList()
for i in model.cities:
    for j in model.cities:
        model.symmetric.add(
            model.x[i,j] == model.x[j,i]
        )

def obj_rule(model):
    # total = 0
    # for i in model.x:
    #     print(model.d[i])
    #     total += model.d[i]*model.x[i]
    # return summation(model.d,model.x)
    return sum(D[i,j]*model.x[i,j] for i in model.cities for j in model.cities)

print("done")

model.obj = Objective(
    rule=obj_rule,
    sense=minimize,
    doc='minimise(distance)'
)

solver = SolverFactory('cplex')
# solver.display.set(2)


results = solver.solve(model, tee = True)

results.write()
