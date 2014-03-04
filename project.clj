(defproject bi-mcmc "0.0.1-SNAPSHOT"
  :description "Bayesian Inference with Markov Chain Monte Carlo"
  :url "http://thegeez.net"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 ;; warning this is a local copy with a bug fix
                 ;; [incanter "1.5.5-SNAPSHOT"]
                 [incanter/incanter-charts "1.5.5-SNAPSHOT"]
                 [incanter/incanter-core "1.5.5-SNAPSHOT"]

                 [org.clojure/data.csv "0.1.2"]]

  :jvm-opts [ "-Xms1500M" "-Xmx1500M"]
  :profiles {:uberjar {:main net.thegeez.bi-mcmc.main
                       :aot [net.thegeez.bi-mcmc.main]}})
