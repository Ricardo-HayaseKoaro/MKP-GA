package problems;

import solutions.Solution;

/**
 * The Evaluator interface gives to a problem the required functionality to
 * obtain a mapping of a solution (n-dimensional array of elements of generic
 * type E (domain)) to a Double (image). It is a useful representation of an
 * objective function for an optimization problem.
 * 
 * @author ccavellucci, fusberti
 * @param <E>
 */
public interface Evaluator<E> {

	/**
	 * Gives the size of the problem domain. Typically this is the number of
	 * decision variables of an optimization problem.
	 * 
	 * @return the size of the problem domain.
	 */
	public abstract Integer getDomainSize();

	/**
	 * The evaluating function is responsible for returning the mapping value of
	 * a solution.
	 * 
	 * @param sol
	 *            the solution under evaluation.
	 * @return the evaluation of a solution.
	 */
	public abstract Double evaluate(Solution<E> sol);

	/**
	 * Evaluates the cost variation of inserting an element into a solution
	 * according to an objective function.
	 * 
	 * @param elem
	 *            the element under consideration for insertion.
	 * @param sol
	 *            the solution for which the element insertion is being
	 *            evaluated.
	 * @return the cost variation resulting from the element insertion into the
	 *         solution.
	 */
	public abstract Double evaluateInsertionCost(E elem, Solution<E> sol);

	/**
	 * Evaluates the cost variation of removing an element into a solution
	 * according to an objective function.
	 * 
	 * @param elem
	 *            the element under consideration for removal.
	 * @param sol
	 *            the solution for which the element insertion is being
	 *            evaluated.
	 * @return the cost variation resulting from the element removal pf the the
	 *         solution.
	 */
	public abstract Double evaluateRemovalCost(E elem, Solution<E> sol);

	/**
	 * Evaluates the cost variation of exchanging candidates, one being
	 * considered to enter the solution (elemIn) and the other being considered
	 * for removal (elemOut).
	 * 
	 * @param elemIn
	 *            the element under consideration for insertion.
	 * @param elemOut
	 *            the element under consideration for removal.
	 * @param sol
	 *            the solution for which the elements exchange is being
	 *            evaluated.
	 * @return the cost variation resulting from the elements exchange.
	 */
	public abstract Double evaluateExchangeCost(E elemIn, E elemOut, Solution<E> sol);

}
