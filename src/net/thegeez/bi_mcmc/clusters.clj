(ns net.thegeez.bi-mcmc.clusters
  (:require [incanter.charts :as charts]
            [incanter.core :as inc-core]
            [incanter.stats :as stats]
            [net.thegeez.bi-mcmc.mcmc :as mcmc]
            [net.thegeez.bi-mcmc.net :as net]
            [net.thegeez.bi-mcmc.distributions :as dists]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.pprint :as pprint]))

;; chapter 3 from
;; http://nbviewer.ipython.org/github/CamDavidsonPilon/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/blob/master/Chapter3_MCMC/IntroMCMC.ipynb

;; from matkelcey
;; from pymc import *
;; from pymc.Matplot import plot
;; #data = [normal(100, 20) for _i in xrange(1000)]
;; data = map(float, open('data', 'r').readlines())
;; >>> min(data)
;; 40.3058046075
;; >>> max(data)
;; 172.745689773
;; mean = Uniform('mean', lower=min(data), upper=max(data))
;; precision = Uniform('precision', lower=0.0001, upper=1.0)
;; process = Normal('process', mu=mean, tau=precision, value=data, observed=True)

;; # do posterior sampling
;; m = pymc.MCMC([mean, precision, process],verbose=10)
;; #m.sample(iter=500, verbose=10)
;; m.sample(iter=10, verbose=10)
;; print(m.stats())

;; >>> m.stats()
;; {'precision': {'95% HPD interval': array([ 0.00307265,  0.03616355]), 'n': 500, 'quantiles': {2.5: 0.0030726520461263188, 25: 0.0056304478695005343, 50: 0.0056304478695005343, 75: 0.0056304478695005343, 97.5: 0.036163554933286587}, 'standard deviation': 0.024166450376244385, 'mc error': 0.0022530354801112954, 'mean': 0.012143625932209764}, 'mean': {'95% HPD interval': array([ 100.00015441,  117.23326122]), 'n': 500, 'quantiles': {2.5: 100.00015440990551, 25: 100.00015440990551, 50: 102.89806412433191, 75: 102.89806412433191, 97.5: 130.95207326510638}, 'standard deviation': 8.3417552405533666, 'mc error': 0.82388925092850029, 'mean': 103.90355067943328}}
;; >>> 


(defonce data (let [file-path "data/matkelcey_simple.csv"]
                (doall
                 (map (comp read-string first) (csv/read-csv (slurp (io/resource file-path)))))))

(def structure {:mean {:type :stochastic
                       :parents []
                       :children [:process]
                       :dist-fn (constantly
                              (dists/uniform-distribution 40.304 172.745)
                              #_(dists/normal-distribution 80.0 20.0)
                              )
                       :n 1
                       }
                :std-dev {:type :stochastic
                          :n 1
                          :parents []
                          :children [:process]
                          :dist-fn (constantly
                                 (dists/uniform-distribution 1.0 50.0)
                                 #_(dists/uniform-distribution 19.9 20.1))}
                :process {:type :stochastic
                          :parents [:mean :std-dev]
                          :children []
                          :dist-fn (fn [{:keys [mean std-dev] :as s}]
                                     (dists/normal-distribution mean std-dev))
                                  
                          :n (count data)
                          :observed data}
                })

(def structure-tau {:mean {:type :stochastic
                           :parents []
                           :children [:process]
                           :dist-fn (constantly
                                  #_(dists/uniform-distribution 99.00 101.00)
                                  (dists/uniform-distribution 40.304 172.745)
                                  #_(dists/normal-distribution 80.0 20.0)
                                  )
                           :n 1
                           }
                    :tau {:type :stochastic
                                :n 1
                                :parents []
                                :children [:std-dev]
                          :dist-fn (constantly
                                 (dists/uniform-distribution 0.0001 1.0)
                                 #_(dists/uniform-distribution 0.0024 0.0026))}
                    :std-dev {:type :deterministic
                              :n 1
                              :parents [:tau]
                              :children [:process]
                              :value-fn (fn [{:keys [tau]}]
                                          (/ 1.0 (Math/sqrt tau)))}
                    :process {:type :stochastic
                              :parents [:mean :std-dev]
                              :children []
                              :dist-fn (fn [{:keys [mean std-dev] :as s}]
                                      (dists/normal-distribution mean std-dev))
                              :n (count data)
                              :observed data}
                    })

(defn go []
  (let [model (mcmc/mcmc structure #_-tau)]
    (println "model" model)
    (let [model (mcmc/sample model :iter 5000)]
      (mcmc/pr-model model)
      (def res model)
      (pprint/pprint (mcmc/stats res)))
    ))

(defn go-tau []
  (let [model (mcmc/mcmc structure-tau)]
    (println "model" model)
    (let [model (mcmc/sample model :iter 5000)]
      (mcmc/pr-model model)
      (def res model))
    ))

(comment
  (let [model res
        mean-vals (->> model
                       :samples
                       (map (comp :value :mean)))
        std-vals (->> model
                       :samples
                       (map (comp :value :precision)))]
    (doto (charts/xy-plot (range 500) mean-vals
                          :y-label "mean")
      (charts/add-points (range 500) std-vals
                         :series-label "std")
          #_(charts/add-function lin -3 3)
          ;; show biggest gains when p is low
          #_(charts/add-function (fn [x]
                          (- (pax_p x)
                             (lin x))) 0 1)
          #_(charts/add-points [from to] [p-from p-to])
          #_(charts/add-points [from to] [logp-from logp-to])
          inc-core/view))

  (let [model res
        mean-vals (->> model
                       :samples
                       (map (comp :value :mean)))]
    (-> (charts/histogram mean-vals
                   :title "mean"
                   :nbins 30
                   :density true
                   )
        #_(add-points [-10 20] [0.01 0.01])
        inc-core/view))

  (let [model res
        mean-vals (->> model
                       :samples
                       (map (comp :value :mean)))
        std-vals (->> model
                       :samples
                       (map (comp :value :precision)))
        ]
    (stats/quantile mean-vals :probs [0.025 0.975]))
  (let [model res
        std-vals (->> model
                       :samples
                       (map (comp :value :precision)))
        ]
    (stats/quantile std-vals :probs [0.025 0.975]))

  (doto (charts/function-plot (fn [s]
                                (/ 1.0 (Math/sqrt s))) 0.001 50.0)
    inc-core/view)

  (doto (charts/function-plot (fn [t]
                                (/ 1.0 (* t t))) 0.001 1.0)
    inc-core/view)
  )

(defonce data-two (let [file-path "data/matkelcey_two.csv"]
                    (doall
                     (map (comp read-string first) (csv/read-csv (slurp (io/resource file-path)))))))


(def structure-two (let [n (count data-two)
                         min (apply min data-two)
                         max (apply max data-two)]
                     {:theta {:type :stochastic
                              :dist-fn (constantly
                                       #_(dists/uniform-distribution 0.31 0.36)
                                       (dists/uniform-distribution 0.0 1.0))
                              :n 1
                              }
                      :bern {:type :stochastic
                             :parents [:theta]
                             :dist-fn (fn [{theta :theta}]
                                        (dists/binomial-distribution 1 theta))
                             :n (count data-two)
                             :step-scale 0.152003093 ;; log(1.0 - 0.1) / log 0.5
                             }
                      :mean1 {:type :stochastic
                              :dist-fn (constantly
                                        (dists/uniform-distribution min max))
                              :n 1}
                      :mean2 {:type :stochastic
                              :dist-fn (constantly
                                        (dists/uniform-distribution min max))
                              :n 1}
                      :std-dev {:type :stochastic
                                :dist-fn (constantly
                                          #_(dists/uniform-distribution 18.0 22.0)
                                          (dists/uniform-distribution 0.0 50.0))
                                :n 1
                                }

                      :mean {:type :deterministic
                             :parents [:bern :mean1 :mean2]
                             :value-fn (fn [{:keys [bern mean1 mean2]}]
                                         ;; (if (= bern 1) mean1 mean2)
                                         (+ (* bern mean1)
                                            (* (- 1 bern) mean2)))
                             :n n}

                      :process {:type :stochastic
                                :parents [:mean :std-dev]
                                :dist-fn (fn [{:keys [mean std-dev]}]
                                           (dists/normal-distribution mean std-dev))
                                :n n
                                :observed data-two}}))

(defn chart-go-two [res]
  (let [model res
          mean1-vals (->> model
                         :samples
                         (map (comp :value :mean1)))
          mean2-vals (->> model
                         :samples
                         (map (comp :value :mean2)))
          std-vals (->> model
                        :samples
                        (map (comp :value :std-dev)))
          theta-vals (->> model
                          :samples
                          (map (comp :value :theta)))] 
      (doto (charts/xy-plot (range) mean1-vals
                            :series-label "mean1"
                            :legend true)
        (charts/add-points (range) mean2-vals
                           :series-label "mean2")
        (charts/add-points (range) std-vals
                           :series-label "std")
        inc-core/view)
      (doto (charts/xy-plot (range) theta-vals
                            :series-label "theta"
                            :legend true)
        inc-core/view)))

(defn go-two []
  (let [model (mcmc/mcmc structure-two)]
    (println "model" model)
    (let [model (time (mcmc/sample model
                                   :iter 100000 :burn 0 #_10000 :thin 1000
                                   ;; :iter 100 :burn 10 :thin 5
                                   ))]
      (mcmc/pr-model model)
      (def res model)
      (pprint/pprint (mcmc/stats res))
      (chart-go-two res)))
    )

(comment
  (chart-go-two res)
  )
;; {:theta {:std-dev 0.013403034589704833, :mean 0.712831870202466, :quantiles {0.975 0.7333746636976199, 0.75 0.7210089328056395, 0.5 0.7143299420175615, 0.25 0.70497601670011, 0.025 0.68798796011754}, :mc-error 7.232812366720091E-4, :n 10000, :accepted 185, :rejected 9814}, :bern {:n 10000, :accepted 0, :rejected 9999}, :mean1 {:std-dev 2.9354356240079986, :mean 132.50650257579827, :quantiles {0.975 135.37761623300153, 0.75 133.3180231268265, 0.5 132.31519770405492, 0.25 131.72547871936558, 0.025 129.23042183191617}, :mc-error 0.12904597314788857, :n 10000, :accepted 106, :rejected 9893}, :mean2 {:std-dev 2.9444003351822747, :mean 136.96123383660253, :quantiles {0.975 141.73971252295715, 0.75 138.8813821712648, 0.5 137.1255005968757, 0.25 134.9752034278358, 0.025 132.8474174171206}, :mc-error 0.17057574621581284, :n 10000, :accepted 251, :rejected 9748}, :std-dev {:std-dev 0.9298255428806178, :mean 49.40049384311787, :quantiles {0.975 49.990456566240994, 0.75 49.79866575538405, 0.5 49.53601803303591, 0.25 49.21113425754308, 0.025 48.18724144586106}, :mc-error 0.050223705760840245, :n 10000, :accepted 209, :rejected 9790}}
