package nature;

public class PairIndex extends PairOwn implements Comparable {

    public int first;
    public int second;

    public PairIndex(double dist, boolean isEdge, int first, int second) {
        super(dist, isEdge);
        this.first = first;
        this.second = second;
    }


}

