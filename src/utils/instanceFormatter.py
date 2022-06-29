from gurobipy import GRB, Model

def solve(knapsack_size, weights, values):
    model = Model()
    model.modelSense = GRB.MAXIMIZE
    model.setParam('OutputFlag', False)

    n_items = len(values) # j
    n_knapsacks = len(knapsack_size) # i

    # x[j] = 1 if j is assigned
    variables = []
    for j in range(n_items):
        # variables.append(model.addVar(vtype=GRB.BINARY))
        variables.append(model.addVar(vtype=GRB.CONTINUOUS, ub = 1, lb = 0))

    model.update()

    # A_ij * x_j <= b[i] for all i
    model.addConstrs(sum(weights[i][j] * variables[j] for j in range(n_items)) <= knapsack_size[i] for i in range(n_knapsacks))

    # max: p_j * x_j
    model.setObjective(sum(variables[j] * values[j] for j in range(n_items)))

    model.optimize()

    shadow_price = model.getAttr(GRB.Attr.Pi)
    return shadow_price

def format_mkpnap1():
    with open('/Users/ricardohk/Dev/GA-Framework/instances/mknap1.txt') as f:
        instances = []
        for line in f:
            if line.strip():
                instances.append(line.split())

        n_tests = int(instances[0][0])
        count = 1
        for n in range(n_tests):
            first_line = instances[count]
            n_vars = int(first_line[0])
            n_constrainst = int(first_line[1])
            optimal = first_line[2]
            count += 1

            p = []
            idx = 0
            size = len(instances[count])
            for j in range(n_vars):
                if (j == size):
                    count += 1
                    idx = 0
                    size += len(instances[count])
                p.append(float(instances[count][idx]))
                idx += 1
            count += 1

            r = []
            for i in range(n_constrainst):
                aux = []
                idx = 0
                size = len(instances[count])
                for j in range(n_vars):
                    if (j == size):
                        count += 1
                        idx = 0
                        size += len(instances[count])
                    aux.append(int(instances[count][idx]))
                    idx += 1
                r.append(aux)
                count += 1
            
            b = []
            for i in range(n_constrainst):
                b.append(int(instances[count][i]))

            count += 1

            shadow_price = solve(b, r, p)
            file_name = '/Users/ricardohk/Dev/GA-Framework/instances/mkp%s.txt' % str(n)
            with open(file_name, 'w') as v:
                for c in first_line:
                    v.write(str(c)+' ')
                v.write('\n')
                for c in p:
                    v.write(str(c)+' ')
                v.write('\n')
                for c in r:
                    for c2 in c:
                        v.write(str(c2)+' ')
                    v.write('\n')
                for c in b:
                    v.write(str(c)+' ')
                v.write('\n')
                for c in shadow_price:
                    v.write(str(c)+' ')
                v.write('\n')

def format_mkpnap2():
    with open('/Users/ricardohk/Dev/GA-Framework/instances/mknap2.txt') as f:
        instances = []
        for line in f:
            if("".join(line.split()).isdigit()):
                instances.append(line.split())

        count = 0
        for n in range(48):
            first_line = instances[count]
            n_vars = int(first_line[1])
            n_constrainst = int(first_line[0])
            count += 1

            p = []
            idx = 0
            size = len(instances[count])
            for j in range(n_vars):
                if (j == size):
                    count += 1
                    idx = 0
                    size += len(instances[count])
                p.append(float(instances[count][idx]))
                idx += 1
            count += 1

            b = []
            idx = 0
            size = len(instances[count])
            for i in range(n_constrainst):
                if (i == size):
                    count += 1
                    idx = 0
                    size += len(instances[count])
                b.append(int(instances[count][idx]))
                idx += 1
            count += 1

            r = []
            for i in range(n_constrainst):
                aux = []
                idx = 0
                size = len(instances[count])
                for j in range(n_vars):
                    if (j == size):
                        count += 1
                        idx = 0
                        size += len(instances[count])
                    aux.append(int(instances[count][idx]))
                    idx += 1
                r.append(aux)
                count += 1
            
        
            optimal = instances[count][0]

            shadow_price = solve(b, r, p)

            first_line.append(optimal)

            first_line[0] = n_vars
            first_line[1] = n_constrainst

            file_name = '/Users/ricardohk/Dev/GA-Framework/instances/mkp%s.txt' % str(n+7)
            with open(file_name, 'w') as v:
                for c in first_line:
                    v.write(str(c)+' ')
                v.write('\n')
                for c in p:
                    v.write(str(c)+' ')
                v.write('\n')
                for c in r:
                    for c2 in c:
                        v.write(str(c2)+' ')
                    v.write('\n')
                for c in b:
                    v.write(str(c)+' ')
                v.write('\n')
                for c in shadow_price:
                    v.write(str(c)+' ')
                v.write('\n')

            count += 1


format_mkpnap1()
format_mkpnap2()