/*----------------------------------------------------------------------
   optimize_histogram.cpp 
    - Optimizes and prints a 1D histogram given a set of data.
 
   Copyright (C) 2015  Trystan Koch

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <algorithm>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_gamma.h>





////////////////////////////////////////////////////////////////////////
/// Prototypes

double LogPosteriorBinWidthProbability(std::vector<double>, int, int,
                                       double, double);
std::vector<double> FindBinErrors(gsl_histogram*, int);
void PrintHistogramAndErrors(gsl_histogram*, std::vector<double>);
void PrintPreviousHeader(std::vector<std::string>);
void PrintHistogramInformation(int);


////////////////////////////////////////////////////////////////////////
/// Implementation Constants

const int MAX_BINS = 200;


////////////////////////////////////////////////////////////////////////
/// Main Function

int main(int argc, char *argv[])
{
  // "header" will contain any comments in the input file, while data
  // is a vector of doubles which we want to make a histogram out of.
  std::vector<std::string> header;
  std::vector<double> data;


  // Read in the data file.
  std::string line;
  while (std::getline(std::cin, line))
  {
    // Assume that any comments are prefaced by a # character
    if (line[0]=='#')
    {
      header.push_back(line);
    }
    else
    {
      data.push_back(std::stod(line));
    }
  }


  // Best to find out what the range of the data now.
  int num_data_points = data.size();
  double data_min = *std::min_element(data.begin(), data.end());
  double data_max = *std::max_element(data.begin(), data.end());


  // Using the method from Knuth (2013), find the best possible
  // number of histogram bins. The best number of bins will have the 
  // greatest posterior probability.
  // Iterates upward through the bin sizes from 1 to MAX_BINS, and sets
  // best_num_bins to the bin number that we will actually use for our
  // data as we output it.
  double current_maximum_log_posterior_probability = 0;
  int best_num_bins = 1;
  for (int i=1;i<=MAX_BINS;i++)
  {
    double new_log_posterior_probability 
        = LogPosteriorBinWidthProbability(data, i, num_data_points, 
                                          data_min, data_max);
    if (current_maximum_log_posterior_probability 
          < new_log_posterior_probability)
    {
      current_maximum_log_posterior_probability 
          = new_log_posterior_probability;
      best_num_bins = i;
    }
  }

  
  // Now that we know the optimal number of bins, use that information
  // to define an optimal histogram.
  gsl_histogram *best_histogram = gsl_histogram_alloc(best_num_bins);
  gsl_histogram_set_ranges_uniform(best_histogram, data_min, data_max);
  

  // Put our data into the histogram.
  for (int i=0; i<num_data_points; i++)
  {
    gsl_histogram_increment(best_histogram, data[i]);
  }


  // From the unnormalized histogram, use equation 46 from Knuth (2013)
  // to find the standard deviation of the probability mass for each 
  // bin. We will output this vector along with the histogram itself.
  // Note that it does not need to be scaled, as the normalization is
  // part of its definition
  std::vector<double> bin_errors 
      = FindBinErrors(best_histogram, num_data_points);

  
  // Normalize the histogram so we find a probability distribution 
  // function.
  double bin_width = (data_max - data_min) / best_num_bins;
  double bin_scale = 1 / (num_data_points * bin_width);
  gsl_histogram_scale(best_histogram, bin_scale);


  // First reprint the previous header
  PrintPreviousHeader(header);

  // Then print some information about the histogram.
  PrintHistogramInformation(best_num_bins);

  // Finally, print out the histogram's data.
  PrintHistogramAndErrors(best_histogram,  bin_errors);

  return 0;
}





////////////////////////////////////////////////////////////////////////
/// Helper Functions


// Knuth (2013) writes a formula for determining the posterior 
// probability -- given a set of data points, an assumption of a 
// piecewise probability function, and the constraint that each piece 
// has equal width -- that the  data points were sampled from the given
// probability distribution.
//
// Because this is exactly the problem of determining which bin width
// to use in a histogram given only the data, we use it here. If we
// choose a number of equal-width bins that give a maximal value for
// this function, it should fit the underlying distribution the best.
//
// Even cursory inspection of these paragraphs is confusing: probability
// should be from 0 to 1, right? I'm not doing the subject justice, but
// it's a convincing paper, and I'm just using the result. It's closer
// to a statistical mechanics argument where you have factorially many
// ways to fit the data points into the bins, but read it for yourself.
// What I know: it gives really nice results for everything I've seen
//
// By far the weirdest place I thought I would encounter the log-gamma
// function: histogram optimization.
double LogPosteriorBinWidthProbability(
    std::vector<double> data,
    int num_bins,
    int num_data_points,
    double data_min,
    double data_max
)
{
  // Because there are multiple long terms to sum, I split the
  // computation into five parts for ease of comprehension.
  double p1, p2, p3, p4, p5;


  // Make a histogram with the number of bins given, equally spaced
  // between the data minimum and the data maximum.
  gsl_histogram* test_histogram = gsl_histogram_alloc(num_bins);
  gsl_histogram_set_ranges_uniform(test_histogram, data_min, 
                                   data_max);


  // Put the data into the histogram.
  for (int i=0; i<num_data_points; i++)
  {
    gsl_histogram_increment(test_histogram, data[i]);
  }


  // Calculates the first four terms of the log of the posterior 
  // probability acording to Equation 31 in Knuth (2013). These terms
  // only depend on the number of data points and the number of bins.
  p1 = num_data_points * std::log(num_bins);
  p2 = gsl_sf_lngamma(num_bins * 0.5);
  p3 = -num_bins * gsl_sf_lngamma(0.5);
  p4 = -gsl_sf_lngamma(num_data_points + (num_bins * 0.5));

  // The fifth term in the posterior probability equation examines how
  // many data points ended up in each bin. This is a sum over the bins.
  p5 = 0;
  for (int i=0; i<num_bins; i++)
  {
    p5 += gsl_sf_lngamma( gsl_histogram_get(test_histogram, i) + 0.5);
  }


  // Clean up
  gsl_histogram_free(test_histogram);


  // Sum the terms in the posterior probability and return that sum.
  // This will be maximal for an optimal bin size.
  return p1 + p2 + p3 + p4 + p5;
}





//
//
//
std::vector<double> FindBinErrors(gsl_histogram* histogram, 
                                  int num_data_points)
{
  unsigned int num_bins = gsl_histogram_bins(histogram);
  double histogram_range = gsl_histogram_max(histogram)
                           - gsl_histogram_min(histogram);

  std::vector<double> bin_errors;
  double data_points_in_bin, term_under_radical, error;

  // This loop iterates over the bins in the optimal histogram we found
  // previously, and uses Knuth's formula to find the standard deviation
  // of the probability mass for each. The program uses this standard 
  // deviation as our error bars.
  for (unsigned int i=0; i<num_bins; i++)
  {
    data_points_in_bin = gsl_histogram_get(histogram, i);

    // There's a square root in the formula, and it's easiest to compute
    // the terms underneath it seperately then taking the root after.
    term_under_radical = data_points_in_bin + 0.5 ;
    term_under_radical *= num_data_points - data_points_in_bin 
                          + 0.5 * ( num_bins - 1 );
    term_under_radical /= num_data_points + 0.5 * num_bins + 1;

    // Finally, apply Knuth's formula.
    error = num_bins / ( histogram_range * ( num_data_points + 0.5 ) );
    error *= std::sqrt(term_under_radical);

    bin_errors.push_back(error);
  }

  return bin_errors;
}





// We'll be printing out four columns of information that will define
// the histogram. Because this task is repetative, I kept it out of the
// main function for readability.
void PrintHistogramAndErrors(gsl_histogram* histogram, 
                             std::vector<double> bin_errors)
{
  // Make sure program puts whole double into the output files.
  // Use scientific notation and constant length for numbers.
  std::cout.precision( std::numeric_limits<double>::max_digits10 );
  std::cout << std::fixed;
  std::cout << std::scientific;

  // Formatting: defines the spacing between data columns
  std::string column_spacing = "    ";

  // For each bin, output the information that Gnuplot uses to create a
  // histogram. This means we need the start and end of the bin, the
  // value of the bin, and the error on the bin.
  double bin_upper_bound, bin_lower_bound;
  for (unsigned int i=0; i<gsl_histogram_bins(histogram); i++)
  {
    gsl_histogram_get_range(histogram, i, &bin_lower_bound, 
                            &bin_upper_bound);
    std::cout << bin_lower_bound << column_spacing;
    std::cout << bin_upper_bound << column_spacing;
    std::cout << gsl_histogram_get(histogram, i) << column_spacing;
    std::cout << bin_errors[i] << std::endl;
  }

  return;
}






// Simple Macro function that just makes the main code easier to read.
void PrintPreviousHeader(std::vector<std::string> header)
{
  for (unsigned int i=0; i<header.size(); i++)
  {
    std::cout << header[i] << std::endl;
  }
}





// Passes information about the histogram we're printing to the output
// file. The meta-information about which run was what is kept so that
// it's easier to see what a file corresponds to.
void PrintHistogramInformation(int best_num_bins)
{
  std::cout << "# BEGIN"                             << std::endl;
  std::cout << "#  Histogram Information"            << std::endl;
  std::cout << "#  "                                 << std::endl;
  std::cout << "#    Optimal Number of Bins: "       << std::endl;
  std::cout << "#      " << best_num_bins            << std::endl;
  std::cout << "#    Columns: "                      << std::endl;
  std::cout << "#      Bin Minimum, Bin Maximum, "
            <<        "Number in Bin, Bin Error"     << std::endl;
  std::cout << "#  "                                 << std::endl;
  std::cout << "# END"                               << std::endl;
}

