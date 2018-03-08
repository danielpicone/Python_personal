# noiteration1.py

from pyomo.environ import *
from pyomo.opt import SolverFactory

# Create a solver
opt = SolverFactory('cplex')

#
# A simple model with binary variables and
# an empty constraint list.
#
model = AbstractModel()
model.n = Param(within=NonNegativeIntegers)
model.J = RangeSet(1,4)
model.c = Param(within=NonNegativeIntegers)
model.x = Var(RangeSet(model.n), within=Binary)
def obj_expression(model):
    return summation(model.c,model.x)
model.o = Objective(rule=obj_expression)
model.constraint = ConstraintList()


model.n = 4
model.c = 4
print(model.c.value)
for i in range(model.n.value):
    print(i)
    model.c[i] = i

# Create a model instance and optimize
instance = model.create_instance()
results = opt.solve(instance)
instance.display()
instance.solutions.load_from(results)


instance.constraint.add( instance.x[2] == 1 )
results = opt.solve(instance)
instance.display()
