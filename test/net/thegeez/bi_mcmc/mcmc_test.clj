(ns net.thegeez.bi-mcmc.mcmc-test
  (:require [clojure.test :refer :all]
            [net.thegeez.bi-mcmc.mcmc :as mcmc]))

(deftest spread-test
  (is (= (mcmc/spread-sample-values {:a {:value 1} :b {:value 2}} 1) [{:a 1 :b 2}]))
  (is (= (mcmc/spread-sample-values {:a {:value 1} :b {:value [2 3]}} 2) [{:a 1 :b 2}{:a 1 :b 3}]))
  (is (= (try (mcmc/spread-sample-values {:a {:value 1} :b {:value [2 3]}} 3)
              (catch Exception e :err)) :err)))
