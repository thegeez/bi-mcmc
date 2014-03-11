(ns net.thegeez.bi-mcmc.sms
    (:require [incanter.charts :as charts]
              [incanter.core :as inc-core]
              [incanter.stats :as stats]
              [net.thegeez.bi-mcmc.mcmc :as mcmc]
              [net.thegeez.bi-mcmc.net :as net]
              [net.thegeez.bi-mcmc.distributions :as dists]
              [clojure.data.csv :as csv]
              [clojure.java.io :as io]
              [clojure.pprint :as pprint]))

;; https://github.com/CamDavidsonPilon/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/blob/master/ExamplesFromChapters/Chapter1/SMS_behaviour.py

(defonce sms-data (let [file-path "data/sms.csv"]
                    (doall
                     (map (comp read-string first) (csv/read-csv (slurp (io/resource file-path)))))))


(def structure (let [alpha (/ 1.0 (stats/mean sms-data))
                     sms-count (count sms-data)]
                 {:lambda1 {:type :stochastic
                            :dist-fn (constantly
                                      (dists/exponential-distribution alpha)
                                      #_(dists/uniform-distribution 15 30)
                                      #_(dists/normal-distribution 19.7 5.0)
                                      #_(dists/uniform-distribution 0 10)
                                      )
                            :n 1}
                  :lambda2 {:type :stochastic
                            :dist-fn (constantly
                                      (dists/exponential-distribution alpha)
                                      #_(dists/uniform-distribution 15 30)
                                      #_(dists/normal-distribution 19.7 5.0)
                                      )
                            :n 1}
                  :tau {:type :stochastic
                        :dist-fn (constantly
                                  (dists/integer-uniform-distribution 0 sms-count))
                        :n 1}
                  :lambda {:type :deterministic
                           :n sms-count
                           :parents [:tau :lambda1 :lambda2]
                           :value-fn (fn [{:keys [tau lambda1 lambda2]}]
                                       (-> []
                                           (into (repeat tau lambda1))
                                           (into (repeat (- sms-count tau) lambda2))))}
                  :observation {:type :stochastic
                                :n sms-count
                                :parents [:lambda]
                                :dist-fn (fn [{:keys [lambda]}]
                                           (dists/poisson-distribution lambda))
                                :observed sms-data}}))
(defn sms-charts [model]
  (let [lambda1-vals (->> model
                          :samples
                          (mapv (comp :value :lambda1)))
        lambda2-vals (->> model
                          :samples
                          (mapv (comp :value :lambda2)))
        tau-vals (->> model
                      :samples
                      (mapv (comp :value :tau)))

        expected-sms-per-day (map (fn [day]
                                    (let [ix (map #(< day %) tau-vals)
                                          l1 (->> (for [[b lb1] (map list ix lambda1-vals)
                                                        :when b]
                                                    lb1)
                                                  (reduce +))
                                          l2 (->> (for [[b lb2] (map list ix lambda2-vals)
                                                        :when (not b)]
                                                    lb2)
                                                  (reduce +))]
                                      (/ (+ l1 l2)
                                         (count tau-vals))))
                                  (range 0 (count sms-data)))
        days (range)]
    (-> (charts/bar-chart days sms-data
                          :x-label "Time (days)"
                          :y-label "Count of text messages received"
                          :series-label "Msgs per day"
                          )
        (charts/add-categories days expected-sms-per-day
                               :series-label "expected sms per day")
        inc-core/view)
    (-> (charts/histogram lambda1-vals
                          :title "lambda1"
                          :nbins 30
                          :density true)
        (charts/set-x-range 15 30)
        (charts/set-y-range 0.0 1.0)
        inc-core/view)
    (-> (charts/histogram lambda2-vals
                          :title "lambda2"
                          :nbins 30
                          :density true)
        (charts/set-x-range 15 30)
        (charts/set-y-range 0.0 1.0)
        inc-core/view)
    (-> (charts/histogram tau-vals
                          :title "tau"
                          :nbins 30
                          :density true)
        (charts/set-x-range 0 70)
        inc-core/view)
))
  


(defn go-sms-data []
    (let [model (mcmc/mcmc structure)]
      (println "model" model)
      (let [model (time (mcmc/sample model
                                     ;; :iter 5000 :burn 2000
                                     :iter 100000 :burn 50000
                                     ))]
        (mcmc/pr-model model)
        (def res model)
        (pprint/pprint (mcmc/stats res))
        (sms-charts model))))

(comment
  (let [vs (repeatedly 1000 #(dists/draw (dists/exponential-distribution (/ 1.0 19.7))))]
    (-> (charts/histogram vs
                          :title "lambda1"
                          :nbins 30
                          :density true)
        (charts/set-x-range 0 100)
        inc-core/view))


  ;; get 3 samples from res (will be sample 50k burn + sample i
  ;; plot histograms based on values so far
  ;; maybe at line of current sample?
  ;; do three steps
  (let [ss (:samples res)]
    (some (fn [s]
            (when (= [:accepted :accepted-by-flip :rejected]
                     [(get-in s [:lambda1 :why])
                      (get-in s [:lambda2 :why])
                      (get-in s [:tau :why])])
              s)) ss))
  (let [[s1 s2 s3] (take 3 (drop 100 (:samples res)))]
    (doall (for [s [s1 s2 s3]]
               (into {}
                     (for [[k v] s
                           :when (contains? v :accepted)]
                       [k (select-keys v [:accepted :accepted-by-flip :rejected :value])]))))
    (let [lambda1-so-far (->> res
                              :samples
                              (take 2000)
                              (mapv (comp :value :lambda1)))
          lambda1-vals (->> res
                            :samples
                            (mapv (comp :value :lambda1)))]
      (-> (charts/histogram lambda1-so-far
                            :title "lambda1 so far"
                            :nbins 30
                            :density true)
          (charts/set-x-range 15 30)
          (charts/set-y-range 0.0 1.0)
          inc-core/view)
      (-> (charts/histogram lambda1-vals
                            :title "lambda1"
                            :nbins 30
                            :density true)
          (charts/set-x-range 15 30)
          (charts/set-y-range 0.0 1.0)
          inc-core/view))
    )
  )
