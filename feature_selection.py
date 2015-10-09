import matplotlib.pyplot as plt, pandas as pd, numpy as np
from sklearn.svm import SVC
from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFECV
# from sklearn.datasets import make_classification


print "reading data"
X = pd.read_csv("combat.csv", index_col=[0,1]).T
y = pd.read_csv("labels.csv", index_col = 0).set_index("gsm_name").sample_class.ix[X.index]

X_subset = X[np.choice(X.columns, 1000)]
print "doing CV"# Create the RFE object and compute a cross-validated score.
svc = SVC(kernel="linear")
# The "accuracy" scoring is proportional to the number of correct
# classifications
rfecv = RFECV(estimator=svc, step=1, cv=StratifiedKFold(y, 2),
              scoring='accuracy')
rfecv.fit(X, y)

print("Optimal number of features : %d" % rfecv.n_features_)

# Plot number of features VS. cross-validation scores
fig = plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (nb of correct classifications)")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
fig.savefig("feature_selection.png")
X.T[rfecv.support_].to_csv("features.csv")