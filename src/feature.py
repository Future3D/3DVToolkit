from abc import ABCMeta,abstractmethod

import numpy as np

class Feature:
    __metaclass__=ABCMeta
    def __init__(self):
        pass

class LocalFeature(Feature):
    __metaclass__=ABCMeta
    def __init__(self,descriptor,rf,point,normal):
        self.descriptor_=descriptor
        self.rf_=rf
        self.point_=point
        self.normal_=normal

    @property
    def descriptor(self):
        return self.descriptor

    @property
    def reference_frame(self):
        return self.rf_

    @property
    def keypoint(self):
        return self.point_
    @property
    def normal(self):
        return self.normal_

class GlobalFeature(Feature):
    __metaclass__ = ABCMeta

    def __init__(self,key,descriptor,pose=None):
        self.key_=key
        self.descriptor_=descriptor
        self.pose_=pose

    @property
    def key(self):
        return self.key_

    @property
    def descriptor(self):
        return self.descriptor_

    @property
    def pose(self):
        return self.pose_

class SHOTFeature(LocalFeature):
    def __init__(self,descriptor,rf,point,normal):
        LocalFeature.__init__(self,descriptor,rf,point,normal)

class MVCNNFeature(GlobalFeature):
    def __init__(self,key,descriptor,pose=None):
        GlobalFeature.__init__(self,key,descriptor,pose)

class BagOfFeature:
    def __init__(self,features=None):
        self.features_=features
        if self.features_ is None:
            self.features_=[]
        self.num_features_=len(self.features_)
    def add_feature(self,feature):
        self.features_.append(feature)
        self.num_features_=len(self.features_)
    def add_features(self,features):
        self.features_.extend(features)
        self.num_features_=len(self.features_)

    def get_feature(self,index):
        if index <0 or index>=self.num_features_:
            raise ValueError('Index %d out of range'%(index))
        return self.features_[index]
    def feature_subset(self,indices):
        if isinstance(indices,np.ndarray)
            indices=indices.tolist()
        if not isinstance(indices,list):
            raise ValueError('Can only index with lists')
        return [self.features_[i] for i in indices]
    @property
    def num_features(self):
        return self.num_features_
    @property
    def descriptors(self):
        return np.array([f.descriptor for f in self.features_])
    @property
    def reference_frame(self):
        return np.array([f.reference_frame for f in self.features_])
    @property
    def keypoints(self):
        return np.array([f.keypoints for f in self.features_])
    @property
    def normals(self):
        return np.array([f.normals for f in self.features_])
