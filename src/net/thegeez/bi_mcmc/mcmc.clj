(ns net.thegeez.bi-mcmc.mcmc
  (:require [incanter.core :as inc-core]
            [incanter.stats :as stats]
            [net.thegeez.bi-mcmc.net :as net]
            [net.thegeez.bi-mcmc.distributions :as dists]
            [clojure.pprint :as pprint]))

(defn mapv2 [f c1 c2]
  (loop [c1 c1
         c2 c2
         res (transient [])]
    (let [c1f (first c1)
          c2f (first c2)]
      (if (and c1f c2f)
        (recur (rest c1)
               (rest c2)
               (conj! res (f c1f c2f)))
        (persistent! res)))))

(defmacro get-inm [m ks]
  `(-> ~m ~@(for [k ks]
            (list 'get k))))

(defmacro assoc-inm [m ks v]
  (if (= (count ks) 2)
    `(assoc ~m ~(first ks) (assoc (get ~m ~(first ks)) ~(second ks) ~v))
    `(assoc-in ~m ~ks ~v)))

(defmacro update-inm [m ks f]
      (if (= (count ks) 2)
        `(assoc-inm ~m ~ks (~f (get-inm ~m ~ks)))
        `(update-in ~m ~ks ~f)))

(defmacro reducem [f init coll]
  `(loop [acc# ~init
          c# ~coll]
     (if-let [c1# (first c#)]
       (recur (~f acc# c1#)
              (rest c#))
       acc#)))

(defn variable-value [variable sample]
  (let [{:keys [n value-fn]} variable]
    (if (< 1 n)
      (loop [i 0
             res (transient [])]
        (if (= i n)
          (persistent! res)
          (or (and (= i 0)
                   (let [first-value (value-fn sample 0)]
                     (when (sequential? first-value)
                       first-value)))
              (recur (inc i)
                     (assoc! res i (value-fn sample i))))))
      #_(mapv (fn [i]
              (value-fn sample i)) (range n))
      (value-fn sample))))

(defn step-value [prev-value prev-step-sd]
  (let [new-value (dists/draw-normal prev-value prev-step-sd)]
    (if (integer? prev-value)
      (Math/round ^double new-value)
      new-value)))

(defn log+ [l r]
  (let [res (+ l r)]
    (if (= res Double/NEGATIVE_INFINITY)
      (reduced res)
     res)))

(defn self-logp [dist-fn sample value]
  (let [logpdf (fn [dists value]
                 (try (Math/log (dists/pdf dists value))
                      (catch Exception e
                        ;; assume the distribution threw an exception
                        ;; on a value outside its reach (at least the
                        ;; case for binomial-distribution with way-off
                        ;; proposals)
                        Double/NEGATIVE_INFINITY
                        )))]
    (if (sequential? value)
      (loop [v value
             i 0
             res 0.0]
        (if-let [v1 (first v)]
          (let [log+ (+ res (logpdf (dist-fn sample i) v1))]
            (if (= log+ Double/NEGATIVE_INFINITY)
              log+
              (recur (rest v)
                    (inc i)
                    log+)))
          res))
      #_(reduce log+
              0.0
              (map-indexed
               (fn [i v]
                 (logpdf (dist-fn sample i) v))
               value)
            )
      (logpdf (dist-fn sample) value))))

(defn logp [structure k sample]
  (let [children (get-inm structure [k :children])
        children-logp
        (loop [c children
               res 0.0]
          (if-let [c1 (first c)]
            (let [log+ (+ res (get-inm sample [c1 :logp]))]
              (if (= log+ Double/NEGATIVE_INFINITY)
                log+
                (recur (rest c)
                       log+)))
            res))
        #_(reduce log+ 0.0 (map #(get-in sample [% :logp]) children))]
    (if (= Double/NEGATIVE_INFINITY children-logp)
      (-> sample
          (assoc-inm [k :self-logp] children-logp)
          (assoc-inm [k :logp] children-logp))
      (let [;;  _ (println "k" k "acc" #_acc "dist" dist "value" value)
            ;;  _ (println "k" k)
            self-logp (if (= (get-inm structure [k :type]) :stochastic)
                        (let [dist-fn (get-inm structure [k :dist-fn])
                              value (get-inm sample [k :value])]
                          (self-logp dist-fn sample value))
                        ;;Deterministic
                        0.0)
            ;;        _ (println "self-logp" self-logp "for" k)
            logp (+ self-logp children-logp)]
        (-> sample
            (assoc-inm [k :self-logp] self-logp)
            (assoc-inm [k :logp] logp))))))

(defn with-dist-fn-n
  "for variables with n = 1 replace dist-fn with a fn that gets the needed values from sample
for variables with n > 1 also replace the dist-fn with a two arity function [sample n] that works for the nth entry"
  [structure]
  (zipmap (keys structure)
          (map
           (fn [{:keys [n dist-fn parents] :as var}]
             (if (= n 1)
               (assoc var :dist-fn (fn wrapped-dist-fn [sample]
                                     (dist-fn
                                      (loop [p parents
                                             res (transient {})]
                                        (if-let [p1 (first p)]
                                          (recur (rest p)
                                                 (assoc! res p1 (get-inm sample [p1 :value])))
                                          (persistent! res)))
                                      #_(zipmap parents
                                                (map #(get-inm sample [% :value]) parents)))))
               (let [single-parents (seq (filter #(= 1 (get-in structure [% :n])) parents))
                     mult-parents (seq (filter #(< 1 (get-in structure [% :n])) parents))]
                 (assoc var :dist-fn (fn wrapped-dist-fn-n [sample n]
                                       (dist-fn
                                        (loop [sp single-parents
                                               mp mult-parents
                                               res (transient {})]
                                          (if-let [sp1 (and sp (first sp))]
                                            (recur (rest sp)
                                                   mp
                                                   (assoc! res sp1 (get-inm sample [sp1 :value])))
                                            (if-let [mp1 (first mp)]
                                              (recur nil
                                                     (rest mp)
                                                     (assoc! res mp1 (get-inm sample [mp1 :value n])))
                                              (persistent! res))))
                                        #_(-> {}
                                                    (into (map (fn [k]
                                                                 [k (get-in sample [k :value])])
                                                               single-parents))
                                                    (into (map (fn [k]
                                                                 [k (get-in sample [k :value n])])
                                                               mult-parents)))))))))
           (vals structure))))

(defn with-value-fn-n
  "for variables with n = 1 replace value-fn with a fn that gets the needed values from sample
for variables with n > 1 also replace the value-fn with a two arity function [sample n] that works for the nth entry"
  [structure]
  (zipmap (keys structure)
          (map
           (fn [{:keys [n value-fn parents] :as var}]
             (if (= n 1)
               (assoc var :value-fn (fn wrapped-value-fn [sample]
                                      (value-fn
                                       (loop [p parents
                                              res (transient {})]
                                         (if-let [p1 (first p)]
                                           (recur (rest p)
                                                  (assoc! res p1 (get-inm sample [p1 :value])))
                                           (persistent! res)))
                                       #_(zipmap parents
                                                        (map #(get-inm sample [% :value]) parents)))))
               (let [single-parents (seq (filter #(= 1 (get-in structure [% :n])) parents))
                     mult-parents (seq (filter #(< 1 (get-in structure [% :n])) parents))]
                 (assoc var :value-fn (fn wrapped-value-fn-n [sample n]
                                        (value-fn
                                           (loop [sp single-parents
                                                  mp mult-parents
                                                  res (transient {})]
                                             (if-let [sp1 (and sp (first sp))]
                                               (recur (rest sp)
                                                      mp
                                                      (assoc! res sp1 (get-inm sample [sp1 :value])))
                                               (if-let [mp1 (first mp)]
                                                 (recur nil
                                                        (rest mp)
                                                        (assoc! res mp1 (get-inm sample [mp1 :value n])))
                                                 (persistent! res))))

                                           #_(-> {}
                                                      (into (map (fn [k]
                                                                   [k (get-in sample [k :value])])
                                                                 single-parents))
                                                      (into (map (fn [k]
                                                                   [k (get-in sample [k :value n])])
                                                                 mult-parents)))))))))
           (vals structure))))

(defn mcmc [structure & {:as args}]
  (let [model (let [s (net/fill-in-structure structure)
                    order (net/order s)
                    update-order (net/update-order s order)
                    with-blankets (net/with-blankets s order)
                    with-dist-fn-n (with-dist-fn-n with-blankets)
                    with-value-fn-n (with-value-fn-n with-dist-fn-n)]
                {:structure with-value-fn-n
                 :order order
                 :update-order (seq update-order)})
        {:keys [order structure]} model
        init (reduce
              (fn [acc k]
                (let [variable (get structure k)]
                  (if (= (:type variable) :deterministic)
                    ;;Deterministic
                    (assoc-in acc [k :value] (variable-value variable acc))
                    ;;Stochastic
                    (let [{:keys [dist-fn observed n parents step-sd step-scale]} variable]
                      (if observed
                        (assoc-in acc [k :value] (vec observed))
                        (if (< 1 n)
                          (let [value (let [first-draw (dists/draw (dist-fn acc 0))]
                                        (if (sequential? first-draw)
                                          first-draw
                                          (mapv (fn [i]
                                                  (dists/draw (dist-fn acc i))) (range n))))]
                            (update-in acc [k] merge {:accepted 0
                                                      :accepted-by-flip 0
                                                      :rejected 0
                                                      :step-sd (or step-sd
                                                                   (if (some zero? value)
                                                                     (repeat (count value)
                                                                             (if (long (first value))
                                                                               1
                                                                               1.0))
                                                                     (map #(Math/abs ^double %) value)))
                                                      :step-scale (or step-scale
                                                                      1.0)
                                                      :value value }))
                          (let [value (dists/draw (dist-fn  (zipmap (keys acc)
                                                                    (map :value (vals acc)))))]
                            (update-in acc [k] merge {:accepted 0
                                                      :accepted-by-flip 0
                                                      :rejected 0
                                                      :step-sd (or step-sd
                                                                   (if (zero? value)
                                                                     (if (long value)
                                                                       1
                                                                       1.0)
                                                                     (Math/abs ^double value)))
                                                      :step-scale (or step-scale
                                                                      1.0)
                                                      :value value}))))))))
              {}
              (apply concat order))
        _ (println "init: " init)
        logps (reduce
               (fn [init k]
                 (logp structure k init))
               init
               (apply concat (reverse order)))]
    (println "logps:" logps)
    (-> model
        (assoc :init init)
        (assoc :first-sample logps))))

(defn single-sample [model last-sample]
  (let [{:keys [structure update-order]} model
      
        ;;        _ (println "update-order: " update-order)
        propose (fn propose [sample k]
                  ;; (println "update k" k)
                (let [variable (get-inm model [:structure k])
                      {:keys [dist-fn n value-blanket logp-blanket]} variable
                      {prev-logp :logp
                       prev-value :value
                       prev-step-sd :step-sd
                       step-scale :step-scale
                       :as prev-sample} (get sample k)
                      ;;_ (println "prev logp" prev-logp " prev-value" prev-value "prev-step-sd" prev-step-sd " step-scale " step-scale)

                      value (if (= k :bern)
                              ;; todo put this in a step method
                              (let [p-jump (- 1.0 (Math/pow 0.5 step-scale))
                                    ;; (/ (Math/log 0.9) (Math/log 0.5)) ;; should be adaptive
                                    ]
                                (mapv (fn [pv]
                                        (if (< (dists/draw-uniform 0.0 1.0) p-jump)
                                          (if (zero? pv) 1 0)
                                          pv)) prev-value))
                              (if (< 1 n)
                                (mapv2 (fn [pv ps]
                                         (step-value pv (* ps step-scale))) prev-value prev-step-sd)
                                (step-value prev-value (* prev-step-sd step-scale))))

                      ;; if self-logp = -Inf then this is useless attempt
                      new-value-sample (assoc-inm sample [k :value] value)
                      self-logp (self-logp dist-fn new-value-sample value)]
                  (if (= self-logp Double/NEGATIVE_INFINITY)
                    ;; if proposed value is impossible no chance
                    ;; of it being better and no need to
                    ;; try to propagate the value
                    (update-inm sample [k :rejected] inc)
                    ;; accept/reject based on propagated value
                    (let [;; update dependent values
                          propagated-values-sample (reducem
                                                    (fn prop-values [s k]
                                                      (let [variable (get structure k)]
                                                        (assoc-inm s [k :value]
                                                                   (variable-value variable s))))
                                                    new-value-sample
                                                    value-blanket)
                          ;; update impacted logp's
                          new-sample (reducem
                                      (fn new-sample [s k]
                                        (logp structure k s))
                                      propagated-values-sample
                                      logp-blanket)
                          logp (get-inm new-sample [k :logp])
                          log-ratio (- logp prev-logp)
                          accepted (and
                                    ;; logp = -Infinity means suggested value is impossible
                                    (not (= logp Double/NEGATIVE_INFINITY))
                                    ;; always acccept new suggestion when it is better
                                    (or (<= 0.0 log-ratio)
                                        ;; if suggestion is less likely
                                        ;; accept by chance of less-likelyness
                                        (< (Math/log ^double (dists/draw-uniform 0.0 1.0))
                                           log-ratio)))]
                      (-> (if accepted
                            (if (pos? log-ratio)
                              (update-inm new-sample [k :accepted] inc)
                              (-> new-sample
                                  (update-inm [k :accepted] inc)
                                  (update-inm [k :accepted-by-flip] inc)))
                            (update-inm sample [k :rejected] inc))
                          ;; add some diagnostics
                          (update-inm [k] (fn [v] (merge v {:proposed value
                                                            :log-ratio log-ratio
                                                            :prev-logp prev-logp
                                                            :prev-value prev-value
                                                            :why (if accepted
                                                                   (if (pos? log-ratio)
                                                                     :accepted
                                                                     :accepted-by-flip)
                                                                   :rejected)})))))))) 
        next #_(loop [uk update-order
                    sample last-sample]
               (if-let [uk1 (first uk)]
                 (recur (rest uk)
                        (propose sample uk1))
                 sample))
        (reducem
         propose
              last-sample
              update-order)]
    next))

(defn replace-coll-val [m]
   (into {} (for [[k v] m]
              [k (into {}
                       (for [[k v] v]
                         [k (if (sequential? v)
                              (str (pr-str (take 3 v)) " ...")
                              v)]))])))
(defn tune [sample-k-back last-sample]
  (println "TUNE!")
  (->> (for [[k v] last-sample]
         (if-let [step-scale (:step-scale v)]
           (let [{prev-accepted :accepted
                  prev-rejected :rejected} (get sample-k-back k)
                  {:keys [accepted rejected]} (get last-sample k)
                  k-accepted (- accepted prev-accepted)
                  k-rejected (- rejected prev-rejected)
                  accept-ratio (double (/ k-accepted
                                          (+ k-accepted k-rejected)))
                  _ (println "tune " k " accept-ratio" accept-ratio "k-acc" k-accepted
                             "k-rej " k-rejected)
                  new-step-scale (cond
                                  ;; might be stuck
                                  (and (zero? accept-ratio)
                                       (< step-scale 0.00001))
                                  1.0

                                  (< accept-ratio 0.001)
                                  (* step-scale 0.1)
                                  (< accept-ratio 0.05)
                                  (* step-scale 0.5)
                                  (< accept-ratio 0.2)
                                  (* step-scale 0.9)
                                  (> accept-ratio 0.95)
                                  (* step-scale 10.0)
                                  (> accept-ratio 0.75)
                                  (* step-scale 2.0)
                                  (> accept-ratio 0.5)
                                  (* step-scale 1.1)
                                  :default step-scale)]
             [k (assoc v :step-scale new-step-scale)])
           [k v]))
       (into {})))

(defn sample [model & {:keys [iter burn thin] :as args
                       :or {burn 0
                            thin 1}}]
  (loop [sample (:first-sample model)
         k-back-sample nil
         samples (object-array (/ (- iter burn) thin)) #_(transient [])
         i 0]
    (when (zero? (mod i 500))
      (println "n: " i))
    (when (zero? (mod i 1000))
      (println "last sample" (replace-coll-val sample)))
    (if (= i iter)
      (assoc model
        :last-sample sample
        :samples (seq samples) #_(persistent! samples))
      (let [next-sample (if (and (< 0 i)
                                 (zero? (mod i 1000)))
                          (tune k-back-sample sample)
                          sample)
            next-sample (single-sample model next-sample)]
        (recur next-sample
               (if (zero? (mod i 1000))
                 next-sample
                 k-back-sample)
               (if (and (<= burn i)
                        (zero? (mod i thin)))
                 (do (aset samples (/ (- i burn) thin) next-sample)
                     samples)
                 samples) #_(conj! samples next-sample)
               (inc i))))))

(defn pr-model [model]
  (let []
    (pprint/pprint (-> model
                       (update-in [:structure] replace-coll-val)
                       (update-in [:init] replace-coll-val)
                       (update-in [:first-sample] replace-coll-val)
                       (update-in [:last-sample] replace-coll-val)
                       (update-in [:samples]
                                  (fn [ss]
                                    (->> (concat (take 10 ss)
                                                 (drop (- (count ss) 10) ss))
                                         (map replace-coll-val))))))))

(defn stats [model]
  (let [vars (->> model
                  :structure
                  (filter (fn [[k v]]
                            (and (= (:type v) :stochastic)
                                 (not (contains? v :observed)))))
                  keys)
        samples (:samples model)]
    (->> (for [v vars]
           [v (let [traces (map (comp :value v) samples)
                    {:keys [accepted rejected]} (get (last samples) v)]
                (merge
                 {:n (count traces)
                  :accepted accepted
                  :rejected rejected}
                 (when-not (< 1 (get-in model [:structure v :n]))
                   {
                    ;; todo the HDP is eq to quantiles only for normal dist?
                    ;; :95pct-HDP-interval (stats/quantile traces :probs [0.025 0.975])
                    :mc-error (let [batches (partition 100 traces)
                                    means (map stats/mean batches)]
                                (/ (stats/sd means)
                                   (Math/sqrt 100)))
                    :quantiles (let [qs [0.025 0.25 0.50 0.75 0.975]]
                                 (zipmap qs
                                         (stats/quantile traces :probs qs)))
                    :mean (stats/mean traces)
                    :std-dev (stats/sd traces)})))])
         (into {}))))
