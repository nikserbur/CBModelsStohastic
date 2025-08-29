#include <Rcpp.h>
#include <cmath>
#include <random>
using namespace Rcpp;

// Вспомогательная функция для вычисления ранга
NumericVector rank_cpp(const NumericVector& x) {
  int n = x.size();
  NumericVector ranks(n);
  std::vector<std::pair<double, int>> values;

  for (int i = 0; i < n; ++i) {
    values.push_back(std::make_pair(x[i], i));
  }

  std::sort(values.begin(), values.end());

  for (int i = 0; i < n; ++i) {
    ranks[values[i].second] = i + 1;
  }

  return ranks;
}

// [[Rcpp::export]]
List cancor_cpp(const NumericMatrix& x, const NumericMatrix& y) {
  int n = x.nrow();
  int p = x.ncol();
  int q = y.ncol();

  if (n <= 1 || p == 0 || q == 0) {
    return List::create(
      _["cor"] = NumericVector::create(0.0),
      _["xcoef"] = NumericMatrix(p, 1),
      _["ycoef"] = NumericMatrix(q, 1)
    );
  }

  // Центрируем данные
  NumericMatrix x_centered(n, p);
  NumericMatrix y_centered(n, q);

  for (int j = 0; j < p; ++j) {
    double mean_x = 0.0;
    for (int i = 0; i < n; ++i) {
      mean_x += x(i, j);
    }
    mean_x /= n;

    for (int i = 0; i < n; ++i) {
      x_centered(i, j) = x(i, j) - mean_x;
    }
  }

  for (int j = 0; j < q; ++j) {
    double mean_y = 0.0;
    for (int i = 0; i < n; ++i) {
      mean_y += y(i, j);
    }
    mean_y /= n;

    for (int i = 0; i < n; ++i) {
      y_centered(i, j) = y(i, j) - mean_y;
    }
  }

  // Вычисляем ковариационные матрицы
  NumericMatrix Sxx(p, p);
  NumericMatrix Syy(q, q);
  NumericMatrix Sxy(p, q);

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += x_centered(k, i) * x_centered(k, j);
      }
      Sxx(i, j) = sum / (n - 1);
    }
  }

  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += y_centered(k, i) * y_centered(k, j);
      }
      Syy(i, j) = sum / (n - 1);
    }
  }

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < q; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += x_centered(k, i) * y_centered(k, j);
      }
      Sxy(i, j) = sum / (n - 1);
    }
  }

  // Проверка на вырожденность
  bool singular = false;
  for (int i = 0; i < p; ++i) {
    if (Sxx(i, i) < 1e-12) {
      singular = true;
      break;
    }
  }
  for (int i = 0; i < q; ++i) {
    if (Syy(i, i) < 1e-12) {
      singular = true;
      break;
    }
  }

  if (singular) {
    return List::create(
      _["cor"] = NumericVector::create(0.0),
      _["xcoef"] = NumericMatrix(p, 1),
      _["ycoef"] = NumericMatrix(q, 1)
    );
  }

  // Фолбэк на R cancor для сложных вычислений
  Function cancor_r("cancor");
  List result = cancor_r(x, y);

  return result;
}

// [[Rcpp::export]]
List RDCSelection_cpp(const DataFrame& data,
                      Nullable<String> target = R_NilValue,
                      double threshold = 0.1,
                      int k_features = 20,
                      int num_permutations = 500) {

  CharacterVector feature_names = data.names();
  int n_features = feature_names.size();
  CharacterVector features;
  NumericVector target_vector;
  bool has_target = false;

  if (target.isNotNull()) {
    String target_name = as<String>(target);
    has_target = true;

    // Находим целевую переменную
    for (int i = 0; i < n_features; ++i) {
      if (feature_names[i] == target_name) {
        target_vector = data[i];
        break;
      }
    }

    // Создаем список признаков без целевой переменной
    for (int i = 0; i < n_features; ++i) {
      if (feature_names[i] != target_name) {
        features.push_back(feature_names[i]);
      }
    }
  } else {
    features = feature_names;
  }

  int n_obs = 0;
  if (features.size() > 0) {
    NumericVector first_feature = data[as<std::string>(features[0])];
    n_obs = first_feature.size();
  }

  NumericVector rdc_scores(features.size());
  CharacterVector rdc_names(features.size());

  // Быстрые C++ генераторы с фиксированным seed
  std::mt19937 gen(123);
  std::uniform_int_distribution<> int_dist;
  std::normal_distribution<> norm_dist(0.0, 1.0);

  // Объявляем R функции один раз для всех циклов
  Function cancor_r("cancor");

  for (int feat_idx = 0; feat_idx < features.size(); ++feat_idx) {
    rdc_names[feat_idx] = features[feat_idx];

    NumericVector x = data[as<std::string>(features[feat_idx])];
    NumericVector y;

    if (has_target) {
      y = target_vector;
    } else {
      // Выбираем случайный другой признак (как в R: sample(features[features != feat], 1))
      CharacterVector other_features;
      for (int i = 0; i < features.size(); ++i) {
        if (i != feat_idx) {
          other_features.push_back(features[i]);
        }
      }
      if (other_features.size() > 0) {
        int_dist = std::uniform_int_distribution<>(0, other_features.size() - 1);
        int random_idx = int_dist(gen);
        y = data[as<std::string>(other_features[random_idx])];
      } else {
        y = x; // Fallback если только один признак
      }
    }

    // Обработка NA/NaN
    NumericVector x_clean(n_obs);
    NumericVector y_clean(n_obs);

    for (int i = 0; i < n_obs; ++i) {
      if (std::isnan(x[i]) || std::isinf(x[i])) {
        x_clean[i] = 0.0;
      } else {
        x_clean[i] = x[i];
      }

      if (std::isnan(y[i]) || std::isinf(y[i])) {
        y_clean[i] = 0.0;
      } else {
        y_clean[i] = y[i];
      }
    }

    // Вычисляем ранги быстро в C++
    NumericVector x_rank = rank_cpp(x_clean);
    NumericVector y_rank = rank_cpp(y_clean);

    int k = std::min(k_features, n_obs - 1);

    // Создаем случайные матрицы быстро в C++
    NumericMatrix x_matrix(n_obs, k);
    NumericMatrix y_matrix(n_obs, k);

    double max_x_rank = *std::max_element(x_rank.begin(), x_rank.end());
    double max_y_rank = *std::max_element(y_rank.begin(), y_rank.end());

    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < n_obs; ++j) {
        double rand_x = norm_dist(gen);
        double rand_y = norm_dist(gen);

        x_matrix(j, i) = std::sin(rand_x * x_rank[j] / max_x_rank);
        y_matrix(j, i) = std::cos(rand_y * y_rank[j] / max_y_rank);
      }
    }

    // Вычисляем каноническую корреляцию
    double canonical_corr = 0.0;
    try {
      List cca_result = cancor_r(x_matrix, y_matrix);
      NumericVector correlations = cca_result["cor"];
      if (correlations.size() > 0) {
        canonical_corr = *std::max_element(correlations.begin(), correlations.end());
      }
    } catch (...) {
      canonical_corr = 0.0;
    }

    if (num_permutations > 0) {
      NumericVector null_scores(num_permutations);

      for (int perm = 0; perm < num_permutations; ++perm) {
        // Быстро перемешиваем y в C++
        NumericVector y_perm = clone(y_clean);
        std::shuffle(y_perm.begin(), y_perm.end(), gen);

        NumericVector y_rank_perm = rank_cpp(y_perm);
        double max_y_rank_perm = *std::max_element(y_rank_perm.begin(), y_rank_perm.end());

        NumericMatrix y_matrix_perm(n_obs, k);
        for (int i = 0; i < k; ++i) {
          for (int j = 0; j < n_obs; ++j) {
            double rand_y_perm = norm_dist(gen);
            y_matrix_perm(j, i) = std::cos(rand_y_perm * y_rank_perm[j] / max_y_rank_perm);
          }
        }

        try {
          List cca_perm = cancor_r(x_matrix, y_matrix_perm);
          NumericVector correlations_perm = cca_perm["cor"];
          if (correlations_perm.size() > 0) {
            null_scores[perm] = *std::max_element(correlations_perm.begin(), correlations_perm.end());
          } else {
            null_scores[perm] = 0.0;
          }
        } catch (...) {
          null_scores[perm] = 0.0;
        }
      }

      // Вычисляем p-value
      int count = 0;
      for (int i = 0; i < num_permutations; ++i) {
        if (canonical_corr >= null_scores[i]) {
          count++;
        }
      }
      double null_mean = (double)count / num_permutations;

      double scaled_score = std::abs(null_mean - 0.5) / 0.5;
      rdc_scores[feat_idx] = std::min(1.0, scaled_score);
    } else {
      rdc_scores[feat_idx] = canonical_corr;
    }
  }

  // Присваиваем имена
  rdc_scores.names() = rdc_names;

  // Отбираем признаки по порогу
  CharacterVector selected_features;
  CharacterVector removed_features;

  for (int i = 0; i < features.size(); ++i) {
    if (rdc_scores[i] >= threshold) {
      selected_features.push_back(features[i]);
    } else {
      removed_features.push_back(features[i]);
    }
  }

  // Добавляем целевую переменную если есть
  CharacterVector final_features;
  if (has_target) {
    for (int i = 0; i < selected_features.size(); ++i) {
      final_features.push_back(selected_features[i]);
    }
    final_features.push_back(as<String>(target));
  } else {
    final_features = selected_features;
  }

  // Создаем фильтрованный DataFrame
  List filtered_data_list;
  for (int i = 0; i < final_features.size(); ++i) {
    std::string feat_name = as<std::string>(final_features[i]);
    filtered_data_list.push_back(data[feat_name]);
  }
  filtered_data_list.names() = final_features;
  DataFrame filtered_data(filtered_data_list);

  return List::create(
    _["data"] = filtered_data,
    _["rdc_scores"] = rdc_scores,
    _["removed_features"] = removed_features,
    _["selected_features"] = final_features,
    _["threshold"] = threshold
  );
}
