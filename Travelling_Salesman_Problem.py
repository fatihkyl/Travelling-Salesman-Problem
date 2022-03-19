import random
import turtle
import math
import copy

def student_details():
    student_id = int("18009814")
    student_username = str("Fatih Ihsan Kayali")
    return student_id, student_username

def generate_map(x_range, y_range, locations):
    generated_map = []
    for x in range(locations):
        random_x_points = random.randint(x_range,y_range)
        random_y_points = random.randint(x_range,y_range)
        generated_map.append([random_x_points, random_y_points])
    return generated_map
#a = generate_map(-200,200,8) uncomment this if print_map doesn't work

def print_map(speed, color, thickness, selected_map):
    print("printing map")
    turtle.penup()
    turtle.setpos(selected_map[0])
    turtle.pendown()
    for x in range(len(selected_map)):
        turtle.speed(speed)
        turtle.pencolor(color)
        turtle.pensize(thickness)
        turtle.goto(selected_map[x])
        turtle.setpos(selected_map[0])
        
#print(print_map(1,'red',2,a)) uncomment this if print_map doesn't work

def calculate_distance(starting_x, starting_y, destination_x, destination_y):
    distance = math.hypot(destination_x - starting_x, destination_y - starting_y)
    return distance


def calculate_path(selected_map):
    distance = 0.0
    current_point = selected_map[0]
    for next_point in selected_map[1:]:
        current_distance = calculate_distance(current_point[0], current_point[1], next_point[0], next_point[1])
        print(current_point, 'to', next_point, '=', current_distance)
        distance += current_distance
        current_point = next_point
    return distance

def nearest_neighbour_algorithm(selected_map):
    def dist(a, b):
        return math.sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2))
    if len(selected_map) == 0:
        return []
    current = selected_map[0]
    optermised_map = [current]
    selected_map.remove(current)
    while len(selected_map) > 0:
        next = selected_map[0]
        for point in selected_map:
            if dist(current, point) < dist(current, next):
                next = point      
        optermised_map.append(next)
        selected_map.remove(next)
        current = next
    return optermised_map

def genetic_algorithm(selected_map, population, iterations, mutation_rate, elite_threshold):
    pop = create_population(population, selected_map)
    print("Initial distance: " + str(1 / fitness_function(pop)[0][1]))
    
    for i in range(0, elite_threshold):
        pop = mating_function(pop,mutation_rate, iterations )
    
    print("Final distance: " + str(1 / fitness_function(pop)[0][1]))
    bestRouteIndex = fitness_function(pop)[0][0]
    bestRoute = pop[bestRouteIndex]
    return bestRoute

def create_population(population, selected_map):
    gene_pool = [selected_map for i in range(0, population)]
    return gene_pool

	
def fitness_function(gene_pool):   
    def path_fitness(cities):
        total_dis = calculate_path(cities)
        fitness= 0.0
        if fitness == 0:
            fitness = 1 / float(total_dis)
        return fitness 
    def rankPathes(population):
        fitnessResults = {}
        for i in range(0,len(population)):
            fitnessResults[i] = path_fitness(population[i])     
        return sorted(fitnessResults.items(),  reverse = True)
    return rankPathes(gene_pool)


#----------------------------------
#-------Helper functions-----------
#----------------------------------

def selection(popRanked, eliteSize):
    selectionResults = []
    df = pd.DataFrame(np.array(popRanked), columns=["Index","Fitness"])
    df['cum_sum'] = df.Fitness.cumsum()
    df['cum_perc'] = 100*df.cum_sum/df.Fitness.sum()
    
    for i in range(0, eliteSize):
        selectionResults.append(popRanked[i][0])
    for i in range(0, len(popRanked) - eliteSize):
        pick = 100*random.random()
        for i in range(0, len(popRanked)):
            if pick <= df.iat[i,3]:
                selectionResults.append(popRanked[i][0])
                break
    return selectionResults

def matingPool(population, selectionResults):
    matingpool = []
    for i in range(0, len(selectionResults)):
        index = selectionResults[i]
        matingpool.append(population[index])
    return matingpool

def breedPopulation(matingpool, eliteSize):
    children = []
    length = len(matingpool) - eliteSize
    pool = random.sample(matingpool, len(matingpool))

    for i in range(0,eliteSize):
        children.append(matingpool[i])
    
    for i in range(0, length):
        child = breed(pool[i], pool[len(matingpool)-i-1])
        children.append(child)
    return children

def mutatePopulation(population, mutationRate):
    mutatedPop = []
    
    for ind in range(0, len(population)):
        mutatedInd = mutate(population[ind], mutationRate)
        mutatedPop.append(mutatedInd)
    return mutatedPop

#--------------------------------------
#-------Helper functions end-----------
#--------------------------------------



def mating_function(gene_pool, mutation_rate, elite_threshold):
    popRanked = fitness_function(gene_pool)
    selectionResults = selection(popRanked, elite_threshold)
    matingpool = matingPool(gene_pool, selectionResults)
    children = breedPopulation(matingpool, elite_threshold)
    new_gene_pool = mutatePopulation(children, mutation_rate)
    return new_gene_pool

def breed(parent_1, parent_2):
    child = []
    childP1 = []
    childP2 = []
    geneA = int(random.random() * len(parent_1))
    geneB = int(random.random() * len(parent_1))
    startGene = min(geneA, geneB)
    endGene = max(geneA, geneB)
    for i in range(startGene, endGene):
        childP1.append(parent_1[i])
    childP2 = [item for item in parent_2 if item not in childP1]
    child = childP1 + childP2
    return child

def mutate(child, mutation_rate):
    for swapped in range(len(child)):
        if(random.random() < mutation_rate):
            swapWith = int(random.random() * len(child))
            city1 = child[swapped]
            city2 = child[swapWith]
            child[swapped] = city2
            child[swapWith] = city1
    return child






