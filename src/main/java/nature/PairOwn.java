package nature;

public class PairOwn implements Comparable {

    public double dist;
    public boolean isEdge;

    public PairOwn(double dist, boolean isEdge) {
        this.dist = dist;
        this.isEdge = isEdge;
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

