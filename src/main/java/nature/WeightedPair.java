package nature;

public class WeightedPair implements Comparable {

    public double dist;
    public double alpha;
    public double beta;


    public WeightedPair(double dist, double alpha, double beta) {
        this.dist = dist;
        this.alpha = alpha;
        this.beta = beta;
    }

    public int compareTo(Object o) {
        PairOwn other = (PairOwn) o;
        if (this.dist < other.dist) {
            return -1;
        }
        if (this.dist > other.dist) {
            return 1;
        }
        return 0;
    }
}

