#include <Rcpp.h>
#include <cmath>
#include <random>
#include <algorithm>
using namespace Rcpp;

// FastMonteCarlo.cpp - Rcpp оптимизации для Monte Carlo симуляций
// Ускорение генерации стохастических процессов
// Автор: Бураков Никита Сергеевич

// [[Rcpp::export]]
NumericMatrix FastAlphaStableSimulation(int n_sims, int horizon, double alpha, 
                                       double beta, double gamma, double delta) {
  NumericMatrix result(n_sims, horizon);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> normal(0, 1);
  std::uniform_real_distribution<> uniform(0, 1);
  
  for (int i = 0; i < n_sims; i++) {
    for (int t = 0; t < horizon; t++) {
      
      // Генерация α-stable распределения через формулу Chambers-Mallows-Stuck
      double u = uniform(gen) * M_PI - M_PI/2;
      double w = -log(uniform(gen));
      
      double stable_value;
      if (alpha == 1.0) {
        stable_value = (2/M_PI) * ((M_PI/2 + beta*u) * tan(u) - 
                                  beta * log((M_PI/2 * w * cos(u))/(M_PI/2 + beta*u)));
      } else {
        double t_alpha = tan(M_PI * alpha / 2);
        double beta_alpha = atan(beta * t_alpha) / alpha;
        stable_value = pow(cos(u - alpha * beta_alpha), -1/alpha) * 
                      cos(u - alpha * beta_alpha + alpha * beta_alpha) * 
                      pow(cos(alpha * beta_alpha), 1/alpha) * 
                      pow(w, (1-alpha)/alpha);
      }
      
      result(i, t) = gamma * stable_value + delta;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix FastCompoundPoissonSimulation(int n_sims, int horizon, double lambda, 
                                           double jump_mean, double jump_sd) {
  NumericMatrix result(n_sims, horizon);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::poisson_distribution<> poisson(lambda);
  std::normal_distribution<> normal(jump_mean, jump_sd);
  
  for (int i = 0; i < n_sims; i++) {
    for (int t = 0; t < horizon; t++) {
      int n_jumps = poisson(gen);
      double jump_sum = 0.0;
      
      for (int j = 0; j < n_jumps; j++) {
        jump_sum += normal(gen);
      }
      
      result(i, t) = jump_sum;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix FastIntegratedStochasticSimulation(int n_sims, int horizon, 
                                                 NumericVector initial_rates,
                                                 NumericMatrix transition_matrix,
                                                 NumericVector regime_means,
                                                 NumericVector regime_sds,
                                                 NumericVector shapley_weights,
                                                 NumericVector inflation_deviation,
                                                 NumericVector gdp_gap,
                                                 NumericVector har_volatility,
                                                 double alpha_stable,
                                                 double gamma_stable,
                                                 double sigma_gaussian,
                                                 double lambda_poisson) {
  
  NumericMatrix result(n_sims, horizon);
  IntegerMatrix regime_paths(n_sims, horizon);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> normal(0, 1);
  std::uniform_real_distribution<> uniform(0, 1);
  std::poisson_distribution<> poisson(lambda_poisson);
  
  int n_states = regime_means.size();
  
  for (int i = 0; i < n_sims; i++) {
    double current_rate = initial_rates[i % initial_rates.size()];
    int current_regime = 0;  // Начинаем с режима 0
    
    for (int t = 0; t < horizon; t++) {
      // Обновление режима
      if (t > 0) {
        double cum_prob = 0.0;
        double rand_val = uniform(gen);
        
        for (int k = 0; k < n_states; k++) {
          cum_prob += transition_matrix(current_regime, k);
          if (rand_val <= cum_prob) {
            current_regime = k;
            break;
          }
        }
      }
      
      regime_paths(i, t) = current_regime;
      
      // Расчет базового drift
      double regime_effect = regime_means[current_regime];
      
      // Взвешенные факторы
      double weighted_inflation = inflation_deviation[std::min(t, (int)inflation_deviation.size()-1)] * 
                                 shapley_weights[0];
      double weighted_gdp = gdp_gap[std::min(t, (int)gdp_gap.size()-1)] * 
                           shapley_weights[1];
      double weighted_volatility = har_volatility[std::min(t, (int)har_volatility.size()-1)] * 
                                  shapley_weights[2];
      
      double base_drift = regime_effect + weighted_inflation + weighted_gdp + weighted_volatility;
      
      // Стохастические компоненты
      double vol_multiplier = sqrt(har_volatility[std::min(t, (int)har_volatility.size()-1)]);
      
      // Гауссовский шок
      double gaussian_shock = normal(gen) * sigma_gaussian * vol_multiplier;
      
      // α-stable скачки
      int n_jumps = poisson(gen);
      double stable_jumps = 0.0;
      
      for (int j = 0; j < n_jumps; j++) {
        double u = uniform(gen) * M_PI - M_PI/2;
        double w = -log(uniform(gen));
        
        double stable_value;
        if (alpha_stable == 1.0) {
          stable_value = (2/M_PI) * ((M_PI/2) * tan(u) - 
                                    log((M_PI/2 * w * cos(u))/(M_PI/2)));
        } else {
          double t_alpha = tan(M_PI * alpha_stable / 2);
          stable_value = pow(cos(u), -1/alpha_stable) * 
                        cos(u) * 
                        pow(w, (1-alpha_stable)/alpha_stable);
        }
        
        stable_jumps += gamma_stable * stable_value * vol_multiplier;
      }
      
      // Обновление ставки
      current_rate = current_rate * 0.95 + base_drift * 0.05 + 
                    gaussian_shock + stable_jumps;
      
      // Ограничения
      current_rate = std::max(0.1, std::min(25.0, current_rate));
      
      result(i, t) = current_rate;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix FastHarVolatilityCalculation(NumericVector volatility_series, 
                                          int daily_window, int weekly_window, 
                                          int monthly_window) {
  int n = volatility_series.size();
  NumericMatrix result(n, 4);  // daily, weekly, monthly, quarterly
  
  for (int i = 0; i < n; i++) {
    // Daily component
    if (i >= daily_window - 1) {
      double daily_vol = 0.0;
      for (int j = 0; j < daily_window; j++) {
        daily_vol += volatility_series[i - j];
      }
      result(i, 0) = daily_vol / daily_window;
    }
    
    // Weekly component
    if (i >= weekly_window - 1) {
      double weekly_vol = 0.0;
      for (int j = 0; j < weekly_window; j++) {
        weekly_vol += volatility_series[i - j];
      }
      result(i, 1) = weekly_vol / weekly_window;
    }
    
    // Monthly component
    if (i >= monthly_window - 1) {
      double monthly_vol = 0.0;
      for (int j = 0; j < monthly_window; j++) {
        monthly_vol += volatility_series[i - j];
      }
      result(i, 2) = monthly_vol / monthly_window;
    }
    
    // Quarterly component (66 days)
    if (i >= 66 - 1) {
      double quarterly_vol = 0.0;
      for (int j = 0; j < 66; j++) {
        quarterly_vol += volatility_series[i - j];
      }
      result(i, 3) = quarterly_vol / 66;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix FastForwardBackwardAlgorithm(NumericVector observations,
                                          NumericMatrix transition_matrix,
                                          NumericVector emission_means,
                                          NumericVector emission_sds,
                                          NumericVector initial_probs) {
  int n_obs = observations.size();
  int n_states = emission_means.size();
  
  NumericMatrix forward_probs(n_obs, n_states);
  NumericMatrix backward_probs(n_obs, n_states);
  
  // Forward pass
  for (int i = 0; i < n_states; i++) {
    double emission_prob = R::dnorm(observations[0], emission_means[i], emission_sds[i], false);
    forward_probs(0, i) = initial_probs[i] * emission_prob;
  }
  
  // Нормализация
  double sum_forward = 0.0;
  for (int i = 0; i < n_states; i++) {
    sum_forward += forward_probs(0, i);
  }
  for (int i = 0; i < n_states; i++) {
    forward_probs(0, i) /= sum_forward;
  }
  
  for (int t = 1; t < n_obs; t++) {
    for (int i = 0; i < n_states; i++) {
      double sum_transition = 0.0;
      for (int j = 0; j < n_states; j++) {
        sum_transition += forward_probs(t-1, j) * transition_matrix(j, i);
      }
      double emission_prob = R::dnorm(observations[t], emission_means[i], emission_sds[i], false);
      forward_probs(t, i) = sum_transition * emission_prob;
    }
    
    // Нормализация
    double sum_forward_t = 0.0;
    for (int i = 0; i < n_states; i++) {
      sum_forward_t += forward_probs(t, i);
    }
    for (int i = 0; i < n_states; i++) {
      forward_probs(t, i) /= sum_forward_t;
    }
  }
  
  // Backward pass
  for (int i = 0; i < n_states; i++) {
    backward_probs(n_obs-1, i) = 1.0;
  }
  
  for (int t = n_obs-2; t >= 0; t--) {
    for (int i = 0; i < n_states; i++) {
      double sum_backward = 0.0;
      for (int j = 0; j < n_states; j++) {
        double emission_prob = R::dnorm(observations[t+1], emission_means[j], emission_sds[j], false);
        sum_backward += transition_matrix(i, j) * emission_prob * backward_probs(t+1, j);
      }
      backward_probs(t, i) = sum_backward;
    }
    
    // Нормализация
    double sum_backward_t = 0.0;
    for (int i = 0; i < n_states; i++) {
      sum_backward_t += backward_probs(t, i);
    }
    for (int i = 0; i < n_states; i++) {
      backward_probs(t, i) /= sum_backward_t;
    }
  }
  
  // Возвращаем объединенные вероятности
  NumericMatrix result(n_obs, n_states);
  for (int t = 0; t < n_obs; t++) {
    double sum_gamma = 0.0;
    for (int i = 0; i < n_states; i++) {
      result(t, i) = forward_probs(t, i) * backward_probs(t, i);
      sum_gamma += result(t, i);
    }
    // Нормализация
    for (int i = 0; i < n_states; i++) {
      result(t, i) /= sum_gamma;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
double FastLogLikelihood(NumericVector observations, 
                        NumericMatrix transition_matrix,
                        NumericVector emission_means,
                        NumericVector emission_sds,
                        NumericVector initial_probs) {
  int n_obs = observations.size();
  int n_states = emission_means.size();
  
  NumericMatrix forward_probs(n_obs, n_states);
  
  // Инициализация
  for (int i = 0; i < n_states; i++) {
    double emission_prob = R::dnorm(observations[0], emission_means[i], emission_sds[i], false);
    forward_probs(0, i) = initial_probs[i] * emission_prob;
  }
  
  double log_likelihood = 0.0;
  
  for (int t = 1; t < n_obs; t++) {
    for (int i = 0; i < n_states; i++) {
      double sum_transition = 0.0;
      for (int j = 0; j < n_states; j++) {
        sum_transition += forward_probs(t-1, j) * transition_matrix(j, i);
      }
      double emission_prob = R::dnorm(observations[t], emission_means[i], emission_sds[i], false);
      forward_probs(t, i) = sum_transition * emission_prob;
    }
    
    // Нормализация и накопление log-likelihood
    double sum_forward = 0.0;
    for (int i = 0; i < n_states; i++) {
      sum_forward += forward_probs(t, i);
    }
    
    log_likelihood += log(sum_forward);
    
    for (int i = 0; i < n_states; i++) {
      forward_probs(t, i) /= sum_forward;
    }
  }
  
  return log_likelihood;
} 