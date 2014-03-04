(ns net.thegeez.bi-mcmc.main
  (:require [net.thegeez.bi-mcmc.clusters :as clusters]
            [net.thegeez.bi-mcmc.sms :as sms])
  (:import [javax.swing JOptionPane])
  (:gen-class))


(defn -main [& args]
  (println "Waiting...")
  (JOptionPane/showMessageDialog nil "Ok to continue")
  (println "continue")
  #_(time (clusters/go-two))
  (time (sms/go-sms-data))
  (println "Done."))
