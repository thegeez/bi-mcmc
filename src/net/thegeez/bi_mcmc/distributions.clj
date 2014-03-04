(ns net.thegeez.bi-mcmc.distributions
  "like incanter.distributions but better suited for mcmc usage patterns and structure describtions"
  (:require [incanter.distributions :as inc-dists])
  (:import [cern.jet.random.tdouble.engine DoubleMersenneTwister]
           [cern.jet.random.tdouble Gamma Beta Binomial DoubleUniform Exponential Normal Poisson]))

;; this code would be much nicer if there was a library with static
;; pdf methods that don't require the random seed initialisation
;; cern.jet only has static cdf methods?

#_(set! *warn-on-reflection* true)

(defprotocol Distribution
  (pdf [d v])
  (draw [d]))

(def ^DoubleUniform static-uniform-distribution (DoubleUniform. (DoubleMersenneTwister.)))

(defn draw-uniform [lower upper]
  (-> (doto static-uniform-distribution
        (.setState lower upper))
      .nextDouble))

(deftype UniformDistribution [lower upper]
  Distribution
  (pdf [this v]
       (-> (doto static-uniform-distribution
             (.setState lower upper))
           (.pdf v)))
  ;; TODO should be (draw [this random-source]) when doing multi threaded
  (draw [this]
        (inc-dists/draw (inc-dists/uniform-distribution lower upper))))

(defn uniform-distribution [lower upper]
  (->UniformDistribution lower upper))

(deftype IntegerUniformDistribution [lower upper]
  Distribution
  (pdf [this v]
       (if (<= lower v upper)
         (double (/ 1 (- upper lower)))
         0.0))
  ;; TODO should be (draw [this random-source]) when doing multi threaded
  (draw [this]
        (long (inc-dists/draw (inc-dists/integer-distribution lower upper)))))

(defn integer-uniform-distribution [lower upper]
  (->IntegerUniformDistribution lower upper))

(def ^Normal static-normal-distribution (Normal. 1.0 1.0 (DoubleMersenneTwister.)))

(defn draw-normal [mean std-dev]
  (-> (doto static-normal-distribution
        (.setState mean std-dev))
      .nextDouble))

(deftype NormalDistribution [mean std-dev]
  Distribution
  (pdf [this v]
       (-> (doto static-normal-distribution
             (.setState mean std-dev))
           (.pdf v)))
  (draw [this]
        (inc-dists/draw (inc-dists/normal-distribution mean std-dev))))

(defn normal-distribution [mean std-dev]
  (->NormalDistribution mean std-dev))

(def ^Binomial static-binomial-distribution (Binomial. 1.0 0.5 (DoubleMersenneTwister.)))

(deftype BinomialDistribution [n p]
  Distribution
  (pdf [this v]
       (-> (doto static-binomial-distribution
             (.setNandP n p))
           (.pdf v)))
  (draw [this]
        (inc-dists/draw (inc-dists/binomial-distribution n p))))

(defn binomial-distribution [n p]
  (->BinomialDistribution n p))

(def ^Exponential static-exponential-distribution (Exponential. 1.0 (DoubleMersenneTwister.)))

(deftype ExponentialDistribution [r]
  Distribution
  (pdf [this v]
       (-> (doto static-exponential-distribution
             (.setState r))
           (.pdf v)))
  (draw [this]
        (inc-dists/draw (inc-dists/exponential-distribution r))))

(defn exponential-distribution [r]
  (->ExponentialDistribution r))

(def ^Poisson static-poisson-distribution (Poisson. 1.0 (DoubleMersenneTwister.)))

(deftype PoissonDistribution [mean]
  Distribution
  (pdf [this v]
       (-> (doto static-poisson-distribution
             (.setMean mean))
           (.pdf v)))
  (draw [this]
        (inc-dists/draw (inc-dists/poisson-distribution mean))))

(defn poisson-distribution [mean]
  (->PoissonDistribution mean))
