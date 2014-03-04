(ns net.thegeez.bi-mcmc.net)

(defn fill-in-structure
  "fill in properties of structure, for input each node must specify :parents and :type :prop when applicable"
  [structure]
  (reduce
   (fn [s [k v]]
     (-> s
         (update-in [k :parents] seq)
         (update-in [k :children] (fnil identity #{}))
         ;; for each child relation set the parent relation
         (as-> s
               (reduce (fn [s c]
                         (update-in s [c :children] (fnil conj #{}) k))
                       s (:parents v)))
         (update-in [k :children] seq)))
   structure
   structure))

(defn order
  "orders the node ids from the structure into a vector of sets of node ids. for every node id all its parents are left of it and its children to the right"
  [structure]
  (loop [s structure
         todo (for [[k v] structure
                    :when (empty? (:parents v))]
                k)
         level 0]
    (if (seq todo)
      (recur (reduce (fn [s k]
                       (assoc-in s [k :level] level)) s todo)
             (->> (select-keys structure todo)
                  vals
                  (mapcat :children)
                  (into #{}))
             (inc level))
      (vec (for [i (range level)]
             (set (for [[k v] s
                        :when (= i (:level v))]
                    k)))))))

(defn update-order [structure order]
  (for [level order
        l level
        :when (and (= (get-in structure [l :type]) :stochastic)
                   (not (contains? (get structure l) :observed)))]
    l))

(defn with-blankets
  "for every stochastic add
   :value-blanket, all the deterministic children for which a value change would cascade in order form
   :logp-blanket, all the stochastic children for which a value change would impact their logp from a value cascade, in reverse order form "
  [structure order]
  (let [keep-in-order (fn [ks order]
                        (let [ks (set ks)]
                          (seq (for [level order
                                     l level
                                     :when (contains? ks l)]
                                 l))))]
    (reduce
     (fn [s k]
       (if (not= (get-in s [k :type]) :stochastic)
         s
         (update-in s [k] merge
                    (loop [vb #{}
                           lb #{k}
                           todo (get-in structure [k :children])]
                      (if (empty? todo)
                        {:value-blanket (seq (keep-in-order vb order))
                         :logp-blanket (seq (reverse (keep-in-order lb order)))}
                        (let [det-children (filter (fn [k]
                                                     (= (get-in structure [k :type]) :deterministic)) todo)
                              stoch-children (filter (fn [k]
                                                       (= (get-in structure [k :type]) :stochastic)) todo)]
                          (recur (into vb det-children)
                                 (-> lb
                                     (into det-children)
                                     (into stoch-children))
                                 (mapcat (fn [k]
                                           (get-in structure [k :children])) det-children))))))))
     structure
     (keys structure))))
