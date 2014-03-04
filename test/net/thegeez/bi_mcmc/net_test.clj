(ns net.thegeez.bi-mcmc.net-test
  (:use [clojure.test])
  (:require [net.thegeez.bi-mcmc.net :as net]
            [clojure.pprint :as pprint]))

(def structure {:a {:parents #{}
                    :children #{:c :d :e}
                    :type :stochastic
                    :value-blanket [:d :e]
                    :logp-blanket [:i :h :g :e :d :c :a]}
                :b {:parents #{}
                    :children #{:h}
                    :type :stochastic
                    :value-blanket []
                    :logp-blanket [:h :b]}
                :c {:parents #{:a}
                    :children #{:f}
                    :type :stochastic
                    :value-blanket []
                    :logp-blanket [:f :c]}
                :d {:parents #{:a}
                    :children #{:g}
                    :type :deterministic}
                :e {:parents #{:a}
                    :children #{:i :h}
                    :type :deterministic}
                :f {:parents #{:c}
                    :children #{}
                    :type :stochastic
                    :value-blanket []
                    :logp-blanket [:f]}
                :g {:parents #{:d}
                    :children #{:i}
                    :type :stochastic
                    :value-blanket []
                    :logp-blanket [:i :g]}
                :h {:parents #{:b :e}
                    :children #{}
                    :type :stochastic
                    :value-blanket []
                    :logp-blanket [:h]}
                :i {:parents #{:e :g}
                    :children #{}
                    :type :stochastic
                    :value-blanket []
                    :logp-blanket [:i]}})

(def structure-order [#{:a :b} #{:c :d :e} #{:f :g :h} #{:i}])

(def structure-update-order [:a :b :c :f :g :h :i])

(def structure-blanket-order [[:a [#{:c :d :e} #{:f :g :h} #{:i}]]
                              [:b [#{:c :d :e} #{:f :g :h} #{:i}]]
                              [:c [#{:f :g :h} #{:i}]]
                              [:d [#{:f :g :h} #{:i}]]
                              [:e [#{:f :g :h} #{:i}]]
                              [:f [#{:i}]]
                              [:g [#{:i}]]
                              [:h [#{:i}]]
                              [:i []]])

(deftest structure-test
  (let [minimal-structure (zipmap (keys structure)
                                  (map (fn [n]
                                         (dissoc n :children)) (vals structure)))
        filled-in-structure (net/fill-in-structure minimal-structure)]
    (is (= filled-in-structure structure))
    (let [order (net/order structure)]
      (is (= order structure-order))
      (let [update-order (net/update-order structure order)]
        (is (= update-order structure-update-order)))
      (let [wo-blankets (zipmap (keys filled-in-structure)
                                (map (fn [n] (dissoc n :value-blanket :logp-blanket))
                                     (vals filled-in-structure)))
            with-blankets (net/with-blankets wo-blankets order)]
        (is (= structure with-blankets))))))

