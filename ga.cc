#include <cmath>
#include <ctime>

#include <algorithm>
#include <bitset>
#include <iostream>
#include <functional>
#include <limits>
#include <random>
#include <sstream>
#include <thread>
#include <vector>

#define NO_VERBOSE = 0
#define LITTLE_VERBOSE 1
#define VERY_VERBOSE 2

using namespace std;

static mt19937_64 rnd;
unsigned func_calls = 0;

/*
  Test functions
 */

double
sphere(const vector<double> &x)
{
  double total = 0.0;
  ++func_calls;
  for (size_t i = 0; i < x.size(); ++i)
    total += x[i] * x[i];
  return total;
}

double
beale(const vector <double> &v)
{
  double x = v[0], y = v[1];

  double ans =
    (1.5 - x + x * y) * (1.5 - x + x * y) +
    (2.25 - x + x * y * y) * (2.25 - x + x * y * y) +
    (2.625 - x + x * y * y * y) * (2.625 - x + x * y * y * y);
  return ans;
}

double
ackley(const vector <double> &v)
{
  const double pi = 3.1415926535897932384626434;
  const double e = 2.7182818284590452353602875;
  double x = v[0], y = v[1];
  double ans = -20 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) -
    exp(0.5 * ( cos(2 * pi * x) + cos(2 * pi * y))) + e + 20;
  return ans;
}

double
rosenbrock(const vector <double> &x)
{
  double ans = 0.0;
  for (unsigned i = 0; i < x.size() - 1; ++i)
    ans += 100.0 * ((x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i])) +
      (x[i] - 1) * (x[i] - 1);
  return ans;
}

/*
  Real number genetic algorithm
 */

struct Range {
  double lo;
  double hi;
};

struct Individual {
  vector<double> genes;
  double result;
};

ostream
&operator<<(ostream &os, Individual const &individual) {
  for (auto it = individual.genes.begin(); it != individual.genes.end(); ++it)
    os << *it << ", ";
  os << " -> " << individual.result;
  return os;
}

static void
initIndividual(Individual & individual,
	       const vector< Range >& ranges)
{
    for (unsigned j = 0; j < ranges.size(); ++j) {
      double d;
      uniform_real_distribution<double> dist(ranges[j].lo,
					     ranges[j].hi);
      d = dist(rnd);
      individual.genes.push_back(d);
    }
}

static void
initIndividuals(vector< Individual >& population,
		unsigned population_size,
		const vector< Range >& ranges)
{
 for (unsigned i = 0; i < population_size; ++i) {
    Individual individual;
    initIndividual(individual, ranges);
    population.push_back(individual);
  }
}


void threaded_section(Individual& individual,
		      function<double(const vector<double> &)> func)
{
  // Function to minimize
  individual.result = func(individual.genes);
}

static void
compareIndividuals(function<double(const vector<double> &)> func,
		   vector< Individual >& population,
		   double& target,
		   bool threaded)
{
  vector<thread> threads(population.size());
  if (threaded) {
    for (unsigned j = 0; j < population.size(); ++j)
      threads[j] = thread(threaded_section, ref(population[j]), func);

    for (unsigned j = 0; j < population.size(); ++j) {
      threads[j].join();
      if (population[j].result < target)
	target = population[j].result;
    }
  } else {
    for (unsigned j = 0; j < population.size(); ++j) {
      population[j].result = func(population[j].genes);
      if (population[j].result < target)
	target = population[j].result;
    }
  }

  // Sort by best
  sort(population.begin(), population.end(),
       [&](const Individual& x, const Individual& y)
       {
	 return x.result < y.result;
       });
}

static void
printPopulation(unsigned generation,
		double target,
		vector< Individual >& population)
{
  for (unsigned j = 0; j < population.size(); ++j)
    cout << "Generation: " << generation << " Individual: " << j << " "
	 << population[j] << endl;

  cout << "Best so far: " << target << endl;
  cout << "Best of generation: " << population[0] << endl;
  cout << "**********************" << endl;
}

static void
makeNewGeneration(unsigned num_parents_kept,
		  unsigned num_parents_for_selection,
		  unsigned num_genes,
		  vector< Individual >& population)
{
  vector< Individual > new_population;
  // Keep best parents
  for (unsigned j = 0; j < num_parents_kept; ++j)
    new_population.push_back(population[j]);

  while (new_population.size() < population.size()) {
    // Select two parents randomly
    uniform_int_distribution<unsigned>
      dist1(0, num_parents_for_selection - 1);
    unsigned first = dist1(rnd);
    unsigned second = dist1(rnd);

    // Select crossover point
    uniform_int_distribution<unsigned> dist2(0, num_genes - 1);
    unsigned crossover = dist2(rnd);
    Individual child1, child2;
    for (unsigned j = 0; j < crossover; ++j) {
      child1.genes.push_back(population[first].genes[j]);
      child2.genes.push_back(population[second].genes[j]);
    }
    // Do crossover point
    uniform_real_distribution<double> dist3(0.0, 1.0);
    double beta = dist3(rnd);
    double p1 = population[first].genes[crossover];
    double p2 = population[second].genes[crossover];
    double res = p1 - beta * (p1 - p2);
    child1.genes.push_back(res);
    res = p2 + beta * (p1 - p2);
    child2.genes.push_back(res);
    for (unsigned j = crossover + 1; j < num_genes; ++j) {
      child1.genes.push_back(population[second].genes[j]);
      child2.genes.push_back(population[first].genes[j]);
    }
    new_population.push_back(move(child1));
    new_population.push_back(move(child2));
  }
  population = move(new_population);
}

static void
mutatePopulation(double mutation_rate,
		 const vector< Range >& ranges,
		 vector< Individual >& population)
{
  unsigned num_mutations = (unsigned)
    round(mutation_rate * population.size() * ranges.size());
  uniform_int_distribution<unsigned>
    dist1(1, population.size() - 1);
  uniform_int_distribution<unsigned>
    dist2(0, ranges.size() - 1);
  for (unsigned j = 0; j < num_mutations; ++j) {
    unsigned pop_index = dist1(rnd);
    unsigned gene_index =  dist2(rnd);
    uniform_real_distribution<double> dist(ranges[gene_index].lo,
					   ranges[gene_index].hi);
    population[pop_index].genes[gene_index] = dist(rnd);
  }
}

// Parameters for geneticAlgorithm function
struct GaParm {
  GaParm(unsigned genes, double lo  = -10.0, double hi = 10.0)
  {
    ranges = vector<Range>(genes, {lo, hi});
  }
  unsigned max_iterations = 10;
  double target = 0.0;
  unsigned population_size = 100;
  unsigned num_parents_kept = 4;
  unsigned num_parents_for_selection = 50;
  vector< Range > ranges;
  unsigned seed = 0;
  double mutation_rate = 0.2;
  unsigned verbose = LITTLE_VERBOSE;
  bool threaded = true;
};

Individual
geneticAlgorithmReal(function<double(const vector<double> &)> func,
		     const GaParm& p)
{
  double target = numeric_limits<double>::max();
  vector< Individual > population;
  clock_t t;

  t = clock();
  rnd.seed(p.seed);
  // Initialize chromosomes
  initIndividuals(population, p.population_size, p.ranges);

  // Optimize
  for (unsigned i = 0; i < p.max_iterations && target > p.target; ++i) {
    // Compare
    compareIndividuals(func, population, target, p.threaded);
    // Output
    if (p.verbose == VERY_VERBOSE)
      printPopulation(target, i, population);
    // Create new generation
    makeNewGeneration(p.num_parents_kept, p.num_parents_for_selection,
		      p.ranges.size(), population);
    // Mutation
    mutatePopulation(p.mutation_rate, p.ranges, population);
  }
  t = clock() - t;
  if (p.verbose == LITTLE_VERBOSE)
    cout << "Time GA: " << t << " clicks (" << ( (float) t) / CLOCKS_PER_SEC
	 << " seconds)." << endl;

  return population[0];
}

/*******************

Bitstring Genetic Algorithm - not fully implemented or tested yet.

********************/

const unsigned num_bits = sizeof(unsigned long) * 8;

struct BitIndividual {
  vector< bitset< num_bits > > genes;
  double result;
};

ostream
&operator<<(ostream &os, BitIndividual const &individual) {
  for (auto it = individual.genes.begin(); it != individual.genes.end(); ++it)
    os << *it << ", ";
  os << " -> " << individual.result;
  return os;
}


static void
initBitIndividuals(vector< BitIndividual >& population,
		unsigned population_size,
		unsigned individual_size)
{
 for (unsigned i = 0; i < population_size; ++i) {
    BitIndividual individual;
    for (unsigned j = 0; j < individual_size; ++j) {
      uniform_real_distribution<double> dist(0.0, 1.0);
      if (dist(rnd) < 0.5)
	individual.genes.push_back(0);
      else
	individual.genes.push_back(1);
    }
    population.push_back(individual);
  }
}

static void
compareBitIndividuals(function<double(const vector<double> &)> func,
		      vector< BitIndividual >& population,
		      double& target,
		      const vector< Range >& ranges)
{
  for (unsigned i = 0; i < population.size(); ++i) {
    vector< double > input;
    for (unsigned j = 0; j < ranges.size(); ++j) {
      double d = (double)
	numeric_limits<unsigned long>::max() / population[i].genes[j].to_ulong();
      input.push_back(d * (ranges[j].hi - ranges[j].lo) + ranges[j].lo);
    }
    population[i].result = func(input);
    if (population[i].result < target)
      target = population[i].result;
  }
  // Sort by best
  sort(population.begin(), population.end(),
       [&](const BitIndividual& x, const BitIndividual& y)
       {
	 return x.result < y.result;
       });
}

static void
printBitPopulation(unsigned generation,
		   double target,
		   vector< BitIndividual >& population)
{
  for (unsigned j = 0; j < population.size(); ++j)
    cout << "Generation: " << generation << " Individual: " << j << " "
	 << population[j] << endl;

  cout << "Best so far: " << target << endl;
  cout << "Best of generation: " << population[0] << endl;
  cout << "**********************" << endl;
}


static void
makeNewBitGeneration(unsigned num_parents_kept,
		     unsigned num_parents_for_selection,
		     unsigned num_genes,
		     vector< BitIndividual >& population)
{
  vector< BitIndividual > new_population;
  // Keep best parents
  for (unsigned i = 0; i < num_parents_kept; ++i)
    new_population.push_back(population[i]);

  while (new_population.size() < population.size()) {
    // Select two parents randomly
    uniform_int_distribution<unsigned>
      dist1(0, num_parents_for_selection - 1);
    unsigned first = dist1(rnd);
    unsigned second = dist1(rnd);

    // Select crossover point
    uniform_int_distribution<unsigned> dist2(0, num_genes - 1);
    uniform_int_distribution<unsigned> dist3(0, num_bits);
    unsigned crossover_var = dist2(rnd);
    unsigned crossover_point = dist3(rnd);
    BitIndividual child1, child2;
    for (unsigned i = 0; i < crossover_var; ++i) {
      child1.genes.push_back(population[first].genes[i]);
      child2.genes.push_back(population[second].genes[i]);
    }
    // Crossover gene
    child1.genes.push_back(population[second].genes[crossover_var]);
    child2.genes.push_back(population[first].genes[crossover_var]);
    for (unsigned i = 0; i < crossover_point; ++i) {
      child1.genes[crossover_var][i] = population[first].genes[crossover_var][i];
      child2.genes[crossover_var][i] =
	population[second].genes[crossover_var][i];
    }
    for (unsigned i = crossover_var + 1; i < num_genes; ++i) {
      child1.genes.push_back(population[second].genes[i]);
      child2.genes.push_back(population[first].genes[i]);
    }
    new_population.push_back(move(child1));
    new_population.push_back(move(child2));
  }
  population = move(new_population);
}

void
mutateBitPopulation(double mutation_rate,
		    unsigned num_genes,
		    vector< BitIndividual >& population)
{
  uniform_real_distribution<double>
    dist(0.0, 1.0);
  for (auto& individual : population) {
    for (unsigned i = 0; i < num_genes; ++i) {
      for (unsigned j = 0; j < num_bits; ++j) {
	if (dist(rnd) < mutation_rate)
	  individual.genes[i][j] = !individual.genes[i][j];
      }
    }
  }
}
BitIndividual
geneticAlgorithmBits(function<double(const vector<double> &)> func,
		     const GaParm& p)
{
  double target = numeric_limits<double>::max();
  vector< BitIndividual > population;

  rnd.seed(p.seed);
  // Initialize chromosomes
  initBitIndividuals(population, p.population_size, p.ranges.size());
  // Compare
  compareBitIndividuals(func, population, target, p.ranges);

  for (unsigned i = 0; i < p.max_iterations || target <= p.target; ++i) {
    // Output
    if (p.verbose == VERY_VERBOSE)
      printBitPopulation(target, i, population);
    // Create new generation
    makeNewBitGeneration(p.num_parents_kept, p.num_parents_for_selection,
			 p.ranges.size(), population);
    // Mutation
    mutateBitPopulation(p.mutation_rate, p.ranges.size(), population);
    // Compare
    compareBitIndividuals(func, population, target, p.ranges);
  }
  return population[0];
}

///////////////////////////////////////////////////

Individual
randomOptimize(function<double(vector<double> &)> func,
	       GaParm& p)
{
  Individual best;
  double target = numeric_limits<double>::max();

  unsigned iterations = p.population_size * p.max_iterations;

  for (unsigned i = 0; i < iterations && target > p.target; ++i) {
    Individual individual;
    initIndividual(individual, p.ranges);
    individual.result = func(individual.genes);
    if (individual.result < target) {
      target = individual.result;
      best = individual;
    }
  }
  return best;
}

#define FUNC sphere
#define GENES 20
#define FUNC_STR(s) STR(s)
#define STR(s) #s

int main(int argc, char *argv[])
{
  Individual individual;

  GaParm p(GENES, -10000.0, 10000.0);
  p.max_iterations = 10000;
  p.population_size = 1000;
  p.mutation_rate = 0.001;
  p.verbose = LITTLE_VERBOSE;

  cout << "Genetic Algorithm Result threaded" << endl;
  individual = geneticAlgorithmReal(FUNC, p);
  cout << individual << endl;
  cout << FUNC_STR(FUNC) ": " << FUNC(individual.genes) << endl;
  double ga_ans_t = FUNC(individual.genes);

  cout << "Genetic Algorithm Result unthreaded" << endl;
  p.threaded = false;
  individual = geneticAlgorithmReal(FUNC, p);
  cout << individual << endl;
  cout << FUNC_STR(FUNC) ": " << FUNC(individual.genes) << endl;
  double ga_ans_u = FUNC(individual.genes);

  cout << "Random Optimization Result" << endl;
  individual = randomOptimize(FUNC, p);
  cout << individual << endl;
  cout << FUNC_STR(FUNC) ": " << FUNC(individual.genes) << endl;
  double rand_ans = FUNC(individual.genes);

  cout << "GA threaded performance: " << rand_ans / ga_ans_t << endl;
  cout << "GA unthreaded performance: " << rand_ans / ga_ans_u << endl;

  return 0;
}
