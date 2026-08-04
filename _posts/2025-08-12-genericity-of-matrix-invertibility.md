---
layout: post
title: Genericity of Matrix Invertibility
category: category
---
In topology, a property that holds on a dense open subset of a space is called a [generic property](https://en.wikipedia.org/wiki/Generic_property) in that space. Intuitively, points satisfying a generic property can be found arbitrarily close to *any* point in the space, and sufficiently small perturbations to such points preserve the property. Loosely speaking, the property holds almost everywhere.[^1]

Matrix invertibility is a generic property in the space of square matrices, meaning the set of all invertible matrices is a dense open subset of the set of all square matrices, $\mathbb{R}^{n \times n}$. 

Here is a proof. For openness, notice that by continuity of the determinant function, for any $U \subseteq \mathbb{R}$ open in $\mathbb{R}$ we must have $\text{det}^{-1}(U)$ open in $\mathbb{R}^{n \times n}$. Since $$\mathbb{R}\setminus\{0\}$$ is open in $\mathbb{R}$, it follows that

$$\text{det}^{-1}(\mathbb{R} \setminus \{0\}) = \{A \in \mathbb{R}^{n \times n} : \det(A) \neq 0\}$$

is open in $\mathbb{R}^{n \times n}$. For density, let $A \in \mathbb{R}^{n \times n}$. If $A$ is invertible, take the constant sequence $$\{A_m\}_{m \in \mathbb{N}} = \{A\}_{m \in \mathbb{N}}$$. If $A$ is noninvertible, observe that 

$$p(\alpha) \equiv \det((1-\alpha) A + \alpha I)$$

is a univariate polynomial in $\alpha$. As $p(1) = \det(I) = 1$, $p(\alpha)$ is not identically zero. The only univariate polynomial that has infinitely many roots is the identically zero polynomial. Thus, $p(\alpha)$ has finitely many roots in $\alpha$ and, therefore, $p(0) = \det(A) = 0$ is an isolated root. Now, for $m \in \mathbb{N}$, consider

$$

p(1/m) = \det\left(\left(1 - \frac{1}{m}\right) A + \frac{1}{m} I\right).

$$

Clearly, $p(1/m) \to p(0) = \det(A) = 0$ as $m \to \infty$. Further, since $\alpha = 0$ is an isolated root, there exists $M \in \mathbb{N}$ such that for all $m > M$, we have $p(1/m) \neq 0$. Therefore, if we set 

$$A_m \equiv \left(1 - \frac{1}{m}\right) A + \frac{1}{m} I$$

and consider the sequence $$\{A_m\}_{m \in \mathbb{N}_{>M}}$$, we obtain a sequence of invertible matrices converging to $A$. Ergo the set of invertible matrices is dense in $\mathbb{R}^{n \times n}$. 

#### An alternative, immediate Zariski-topological proof

The proof above works under the standard Euclidean topology. Alternatively, we can study the topological properties of invertible matrices by identifying $\mathbb{R}^{n \times n}$ with $\mathbb R^{n^2}$ and equipping this space with the [Zariski topology](https://en.wikipedia.org/wiki/Zariski_topology), in which closed sets are defined by polynomial equations; that is, by [algebraic varieties](https://en.wikipedia.org/wiki/Algebraic_variety) in $\mathbb R^{n^2}$.

Since the determinant is a polynomial function, the set

$$\{A\in \mathbb{R}^{n \times n}:\det(A) = 0\}$$

of noninvertible matrices is an algebraic variety and hence Zariski closed. Consequently, its complement, which is the set of all invertible matrices in $\mathbb{R}^{n \times n}$,

$$\{A\in \mathbb{R}^{n \times n}:\det(A)\neq 0\},$$

is Zariski open.

Every nonempty Zariski-open subset of $\mathbb R^{n^2}$ is automatically Zariski dense.[^2] Therefore, the set of invertible matrices is Zariski open and dense in $\mathbb{R}^{n \times n}$. Moreover, over the reals, every Zariski-open dense set is also open and dense in the usual Euclidean topology.[^3] Thus, once $\mathbb R^{n\times n}$ is equipped with the Zariski topology, the genericity of matrix invertibility follows immediately from the observation that the set of noninvertible matrices is an algebraic variety.

I find this an interesting, concrete example of how a change of topology can simplify a problem. Once the relation between the Zariski and Euclidean topologies is understood, the genericity of matrix invertibility becomes nearly immediate under the Zariski topology. As a related example, passing to the Zariski topology also yields an elegant proof of the [Cayley–Hamilton theorem](https://en.wikipedia.org/wiki/Cayley%E2%80%93Hamilton_theorem).[^4]


---
{: data-content="links"}

- [Generic property](https://en.wikipedia.org/wiki/Generic_property)
- [Algebraic variety](https://en.wikipedia.org/wiki/Algebraic_variety)
- [Zariski topology](https://en.wikipedia.org/wiki/Zariski_topology)
- [Why Zariski topology?](https://math.stackexchange.com/questions/161884/why-zariski-topology)
- [What is the Zariski topology good/bad for?](https://mathoverflow.net/questions/21502/what-is-the-zariski-topology-good-bad-for)
- [Why is the Zariski topology coarser than standard topology](https://math.stackexchange.com/questions/4093055/why-is-the-zariski-topology-coarser-than-standard-topology)
- [Zariski Open Sets are Dense?](https://math.stackexchange.com/questions/179652/zariski-open-sets-are-dense)
- [Andrew Paul's notes on the Zariski Topology](https://anpaul.onrender.com/myassets/docs/BlogDocs/ZariskiTopology.pdf) 
- [Cayley–Hamilton theorem](https://en.wikipedia.org/wiki/Cayley%E2%80%93Hamilton_theorem)
- [Using Zariski Topology to Prove the Cayley–Hamilton Theorem](https://www.aergus.net/blog/posts/using-zariski-topology-to-prove-the-cayley-hamilton-theorem.html)


---
{: data-content="footnotes"}


[^1]: Genericity can be thought of as the topological counterpart of the measure-theoretic notion of holding almost everywhere. The two notions, however, are not equivalent in general. A property may be generic without holding almost everywhere, or hold almost everywhere without being generic. Nevertheless, the two notions may coincide in particular settings, as they do for matrix invertibility in $\mathbb{R}^{n \times n}$, discussed in this post.

[^2]: See Proposition 6 in [these notes](https://anpaul.onrender.com/myassets/docs/BlogDocs/ZariskiTopology.pdf) by [Andrew Paul](https://anpaul.onrender.com/).

[^3]: The converse is not true. The Zariski topology is coarser than the standard Euclidean topology. See Proposition 4 in Andrew Paul's [notes](https://anpaul.onrender.com/myassets/docs/BlogDocs/ZariskiTopology.pdf). See also [this discussion](https://math.stackexchange.com/questions/4093055/why-is-the-zariski-topology-coarser-than-standard-topology).

[^4]: See [this](https://www.aergus.net/blog/posts/using-zariski-topology-to-prove-the-cayley-hamilton-theorem.html). Expanded details of this proof can also be found in Section 3 of Andrew Paul's [notes](https://anpaul.onrender.com/myassets/docs/BlogDocs/ZariskiTopology.pdf).

