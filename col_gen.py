def create_branches(node_i,constraints):
    branches = [node(),node()]
    branches[0].create_node(node_i.c,node_i.eqmat,node_i.eqb,node_i.D,node_i.d,consets,varsets,
                            A = node_i.A,rhs = node_i.rhs,
                            sub_branch_matrix = node_i.sub_branch_matrix,
                            sub_branch_rhs = node_i.sub_branch_rhs)
    branches[0].add_branch_constraints(constraints[0])
    branches[1].create_node(node_i.c,node_i.eqmat,node_i.eqb,node_i.D,node_i.d,consets,varsets,
                            A = node_i.A,rhs = node_i.rhs,
                            sub_branch_matrix = node_i.sub_branch_matrix,
                            sub_branch_rhs = node_i.sub_branch_rhs)
    branches[1].add_branch_constraints(constraints[1])
    # Add new nodes to the list
    node_list.extend(branches)
    # Delete the old node from the list
    node_list.pop(0)

def solve_node(node):
    node.create_node(model.c,model.eqmat,model.eqb,model.D,model.d,consets,varsets)
    node.column_generation()
    node.get_solution_values()

    for k in node.solution:
            print(node.solution[k])
            for index,point in enumerate(node.solution[k]):
                print("Point is: {0}".format(point))
                print(abs(point-np.round(point)))
                if max_fractional[0] < abs(point-np.round(point)):
                    max_fractional[0] = abs(point-np.round(point))
                    max_fractional[1] = index
                if max_fractional[0]==0.5:
                    break
            if max_fractional[0]==0.5:
                    break

    return node.solution

from __future__ import division
from pyomo.environ import *
from pyomo.core.base.expr import identify_variables
import numpy as np
import logging
logging.getLogger('pyomo.core').setLevel(logging.ERROR)

model = ConcreteModel()


model.n = range(6)

model.x = Var(model.n, within=Binary)

# model.c = np.array([1,2,2,1])
model.c = np.array([1,2,-1,3,-2,1])
# model.a = np.matrix([[-1,-1,0,0],[0,0,-1,-1]])
model.D = np.matrix([[-1,-4,1,0,0,0],[0,0,0,-3,-2,-1]])
# model.b = np.array([-1,-2])
model.d = np.array([-1,-2])
# model.eqmat = np.matrix([1,1,1,1])
model.eqmat = np.matrix([2,-5,1,-2,1,1])
# model.eqb = np.array([3])
model.eqb = np.array([0])
model.set = range(len(model.d))
consets = [[0],[1]]
# varsets = [[0,1],[2,3]]
varsets = [[0,1,2],[3,4,5]]
num_of_sets = len(consets)

def projectBinary(x):
    for row in range(len(x)):
        if x[row] <= 0.5:
            x[row] = 0
        else:
            x[row] = 1
    return x

def projectHyperplane(x,A,b):
    for row in range(A.shape[0]):
        if np.dot(A[row,:],x) - b[row] > 0:
            x = x - ((( (np.dot(np.array(A[row,:]),x))[0] - b[row] )/ ((np.array(np.square(A).sum(axis=1))[row])[0]) ) * np.array(A[row,:].flatten())[0])
    return x

def find_points(A,b):
    """Finds feasible points for a set of linear constraints"""
    points = []
    for i in range(1000):
        counter = 0
        x = np.random.rand(A.shape[1])
        x_old = x+1
        while (np.linalg.norm(x-x_old) > 0.01):
            if counter > 10:
                break
            x_old = x
            x = projectBinary(x)
            x = projectHyperplane(x,A,b)
            counter += 1
        if counter < 10:
            points.append(x)

    return np.unique(np.array(points),axis=0)


def create_sets(consets,varsets,Aineq,bineq):
    num_of_sets = len(consets)
    sets = {}
    lambda_sets = {}
    counter = 0
    for k in range(num_of_sets):
        set_k = {}
        set_k['constraint_rows'] = consets[k]
        set_k['indicies'] = varsets[k]
        set_k['points'] = find_points(Aineq[set_k['constraint_rows'],set_k['indicies']],bineq[set_k['constraint_rows']])
#         set_k['points'] = points1[k]
        set_k['num_points'] = len(set_k['points'])
        for count in range(set_k['num_points']):
            lambda_sets[counter] = (k,count)
            counter += 1
        sets[k] = set_k
    sets['num_sets'] = num_of_sets
    return sets,lambda_sets

points1 = [np.array([[1,1],[0,1]]),np.array([[1,1]])]

sets,lambda_sets = create_sets(consets,varsets,model.D,model.d)

def create_A_and_rhs(eqmat,eqb,sets):
    num_sets = sets['num_sets']
    num_eq = len(eqb)
    big_matrix = np.zeros((num_eq+num_sets,1))
    for k in range(num_sets):
        matrix_of_set_k = np.zeros((num_eq+num_sets,sets[k]['num_points']))
        # Now populate the matrix
        matrix_of_set_k[num_eq+k,:] = 1
        for j in range(num_eq):
            for column,point in enumerate(sets[k]['points']):
                matrix_of_set_k[j,column] = np.dot(eqmat[j,sets[k]['indicies']],point)[0,0]
        big_matrix = np.concatenate((big_matrix,matrix_of_set_k),axis=1)
    big_matrix = np.delete(big_matrix,0,axis=1)

    rhs = np.concatenate((model.eqb,np.ones(num_sets)))
    return big_matrix,rhs


A,rhs = create_A_and_rhs(model.eqmat,model.eqb,sets)

def get_obj_values(c,sets):
    num_sets = sets['num_sets']
    obj_array = np.zeros(1)
    for k in range(num_sets):
        for point in sets[k]['points']:
            obj = np.dot(c[sets[k]['indicies']],point)
            obj_array = np.append(obj_array,obj)
    obj_array = np.delete(obj_array,0)
    return obj_array


def get_new_column_values(k,sets,point,eqmat,eqb,c):
    num_sets = sets['num_sets']
    num_eq = len(eqb)
    obj = np.dot(c[sets[k]['indicies']],point)
    Acol = np.zeros((num_eq+num_sets,1))
    Acol[k+num_eq,0] = 1
    for i in range(num_eq):
        Acol[i,0] = np.dot(eqmat[i,sets[k]['indicies']],point)[0,0]
    return obj,Acol

class node:

#     @classmethod
    def create_node(self,c,eqmat,eqb,D,d,consets,varsets,A = None,rhs = None,sub_branch_matrix = None,sub_branch_rhs = None, incumbent = np.inf):
        self.c = c
        self.D = D
        self.d = d
        self.eqmat = eqmat
        self.eqb = eqb
        self.varsets = varsets
        self.sets,self.lambda_sets = create_sets(consets,varsets,D,d)
        if A is None:
            self.A,self.rhs = create_A_and_rhs(eqmat,eqb,sets)
        else:
            self.A = A
            self.rhs = rhs
        self.obj = get_obj_values(c,self.sets)

        self.sub_branch_matrix = sub_branch_matrix
        self.sub_branch_rhs = sub_branch_rhs

#     @classmethod
    def add_branch_constraints(self,branch_con):
        k = branch_con[0]
        index = branch_con[1]
        value = branch_con[2]
        points = sets[k]['points']
        # find all lambdas in set k
        constraint_matrix = np.ndarray((0,len(lambda_sets)))
        for i in lambda_sets:
            if lambda_sets[i][0]==k:  # Ensures we're in the correct set
                if points[lambda_sets[i][1]][index]!=value:
                    row = np.zeros(len(lambda_sets))
                    row[i] = 1
                    constraint_matrix = np.concatenate((constraint_matrix,[row]), axis = 0)

        constraint_rhs = np.zeros(constraint_matrix.shape[0])
        self.A = np.concatenate((self.A,constraint_matrix),axis=0)
        self.rhs = np.concatenate((self.rhs,constraint_rhs))

        constraint_row = np.zeros((1,len(self.c)))
        constraint_row[0,self.sets[k]['indicies'][index]] = 1
        if self.sub_branch_matrix is not None:
            self.sub_branch_matrix = np.concatenate((self.sub_branch_matrix,constraint_row),axis=0)
            self.sub_branch_rhs = np.concatenate((self.sub_branch_rhs,[value]))
        else:
            self.sub_branch_matrix = np.array((constraint_row))
            self.sub_branch_rhs = np.array([value])

#     @classmethod
    def column_generation(self):
        rmlp = ConcreteModel()
        num_lambda = range(self.A.shape[1])
        rmlp.n = num_lambda
        rmlp.lambda_var = Var(rmlp.n,within=NonNegativeReals)
        rmlp.con = ConstraintList()
        for i in range(A.shape[0]):
            rmlp.con.add(expr = sum( self.A[i,j]*rmlp.lambda_var[j] for j in rmlp.n ) == self.rhs[i])

        rmlp.OBJ = Objective(expr =
                            sum( self.obj[j]*rmlp.lambda_var[j] for j in rmlp.n ))

        rmlp.dual = Suffix(direction=Suffix.IMPORT)
        solver = SolverFactory('cplex',solver_io='nl')
        solver.options['presolve'] = 0
        solver.options['lpdisplay'] = 0
        results = solver.solve(rmlp, tee = False)
        rmlp.solutions.load_from(results)

        convergenceflag = 0
        while convergenceflag == 0:
            convergenceflag = 1
            num_sets = sets['num_sets']
            from pyomo.core import Constraint
            pi = np.zeros(len(self.eqb))
            mu = np.zeros(num_sets)
            for c in rmlp.component_objects(Constraint, active=True):
                cobject = getattr(rmlp, str(c))
                for index in cobject:
                    if index <= len(self.eqb):
                        pi[index-1] = rmlp.dual[cobject[index]]
                    else:
                        mu[index-1-len(pi)] = rmlp.dual[cobject[index]]

            # Now create the subproblem and use the dual values to solve

            for k in range(num_sets):
                num_points = sets[k]['num_points']
                subprob = ConcreteModel()
                subprob.x = Var(range(len(sets[k]['indicies'])),within=Binary)
                subprob.con = ConstraintList()
                for set_num,i in enumerate(sets[k]['constraint_rows']):
                    subprob.con.add(expr =
                                   sum(self.D[i,j]*subprob.x[col] for col,j in enumerate(sets[k]['indicies'])) <= self.d[i]
                                   )

                if self.sub_branch_matrix is not None:
                    subprob.branch_constraints = ConstraintList()
                    for i in range(len(self.sub_branch_rhs)):
                        subprob.branch_constraints.add(expr =
                                                      sum( self.sub_branch_matrix[i,j]*subprob.x[j] for j in range(len(sets[k]['indicies'])))
                                                       == self.sub_branch_rhs[i])

                constants = self.c[sets[k]['indicies']] - np.matmul(pi,np.array(self.eqmat[:,sets[k]['indicies']]))
                subprob.OBJ = Objective(expr =
                                       sum(constants[i]*subprob.x[i] for i in range(len(sets[k]['indicies']))) - mu[k]
                                       )
                solver = SolverFactory('cplex')
                results = solver.solve(subprob, tee = False)
                subprob.solutions.load_from(results)
                if subprob.OBJ() < -10**-10:
                    convergenceflag = 0
                    xval = []
                    for v in subprob.component_objects(Var, active=True):
                        for index in v:
                            xval.append(v[index].value)
                    print(sets[k]['points'])
                    print(np.array(xval))
                    sets[k]['points'] = np.append(sets[k]['points'],[np.array(xval)],axis=0)
                    lambda_sets[np.int(A.shape[1])] = (k,len(sets[k]['points'])-1)
                    newc,newAcol = get_new_column_values(k,sets,xval,model.eqmat,model.eqb,model.c)
                    self.obj = np.append(self.obj,newc)
                    self.A = np.append(self.A,newAcol,axis=1)

            rmlp = ConcreteModel()

            num_lambda = range(self.A.shape[1])
            rmlp.n = num_lambda
            rmlp.lambda_var = Var(rmlp.n,within=NonNegativeReals)
            rmlp.con = ConstraintList()
            for i in range(self.A.shape[0]):
                rmlp.con.add(expr = sum( self.A[i,j]*rmlp.lambda_var[j] for j in rmlp.n ) == self.rhs[i])

            rmlp.OBJ = Objective(expr =
                                sum( self.obj[j]*rmlp.lambda_var[j] for j in rmlp.n ))

            rmlp.dual = Suffix(direction=Suffix.IMPORT)
            solver = SolverFactory('cplex',solver_io='nl')
            solver.options['presolve'] = 0
            solver.options['lpdisplay'] = 0
            results = solver.solve(rmlp, tee = False)
            rmlp.solutions.load_from(results)

        for v in rmlp.component_objects(Var, active=True):
            output = []
            for index in v:
                output.append(v[index].value)
        self.lp_relax_value = rmlp.OBJ()
        self.lp_relaxation = np.array(output)
        self.term_cond = results.solver.termination_condition

    def get_solution_values(self):
        solution = {}
        for index,i in enumerate(self.lp_relaxation):
            if i != 0:
                if self.lambda_sets[index][0] in solution:
                    solution[lambda_sets[index][0]] += i*self.sets[self.lambda_sets[index][0]]['points'][self.lambda_sets[index][1]]
                else:
                    solution[self.lambda_sets[index][0]] = i*self.sets[self.lambda_sets[index][0]]['points'][self.lambda_sets[index][1]]
        self.solution = solution

    def create_branch_constraints(self):
        # Call this to determine if the solution is fractional and if so, create branching sets
        max_fractional = [0,None]
        for k in self.solution:
            for index,point in enumerate(self.solution[k]):
                if max_fractional[0] < abs(point-np.round(point)):
                    max_fractional[0] = abs(point-np.round(point))
                    max_fractional[1] = index
                if max_fractional[0]==0.5:
                    break
            if max_fractional[0]==0.5:
                    break
        if max_fractional[1] is None:
            self.solution['fractional'] = 'n'
            return None
        self.index_branch_var = max_fractional[1]
        num_sets = self.sets['num_sets']
        for k in range(num_sets):
            if self.index_branch_var in self.sets[k]['indicies']:
                return (k,sets[k]['indicies'].index(self.index_branch_var),0),(k,sets[k]['indicies'].index(self.index_branch_var),1)

node_list = []
incumbent_value = np.inf
incumbent_solution = None
best_node = None
node_list.append(node())
node_list[0].create_node(model.c,model.eqmat,model.eqb,model.D,model.d,consets,varsets)
print("Node     Nodes Left       relaxation     incumbent       gap")

node_number = 0
while len(node_list)!=0:
    node_i = node_list[0]
    node_i.column_generation()
    node_i.get_solution_values()
    constraints = node_i.create_branch_constraints()
    if str(node_i.term_cond)=='infeasible':
        # Prune due to infeasibility
        node_list.pop(0)
        status = 'infeasible'
    elif constraints is None:
        if node_i.lp_relax_value < incumbent_value:
            incumbent_value = node_i.lp_relax_value
            incumbent_solution = node_i.lp_relaxation
            best_node = node_i
        # Prune due to optimality
        node_list.pop(0)
        status = 'cutoff'
    else:
        # Create branches
        create_branches(node_i,constraints)
        status = node_i.lp_relax_value
    gap = (node_i.lp_relax_value-incumbent_value)/(incumbent_value+0.0001)
    if np.isnan(gap):
        gap=100
    if incumbent_value==np.inf:
        incumbent_value_print = ''
    else:
        incumbent_value_print = incumbent_value
    nodes_left = len(node_list)
    print_output = [str(node_number),str(nodes_left),str(status),str(incumbent_value_print),str(gap)+'%']
#     print("{:}    {:}     {:}    {:}%".format(*print_output).rjust(10))
    print "".join(word.ljust(14) for word in print_output)
    node_number+=1

if incumbent_solution is None:
    print("The problem is infeasible")
else:
    print("The final solution is {0} with value {1}".format(incumbent_solution,incumbent_value))
    print("The variable values from the OG problem are: {0}".format(best_node.solution))
