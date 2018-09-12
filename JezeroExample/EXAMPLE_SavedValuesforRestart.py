np.save('XconcC.npy', self.domain.quantities['concentration'].centroid_values)
np.save('XelevC.npy', self.domain.quantities['elevation'].centroid_values)
np.save('XxmC.npy', self.domain.quantities['xmomentum'].centroid_values)
np.save('XymC.npy', self.domain.quantities['ymomentum'].centroid_values)
np.save('XstageC.npy', self.domain.quantities['stage'].centroid_values)
np.save('XconcV.npy', self.domain.quantities['concentration'].vertex_values)
np.save('XelevV.npy', self.domain.quantities['elevation'].vertex_values)
np.save('XxmV.npy', self.domain.quantities['xmomentum'].vertex_values)
np.save('XymV.npy', self.domain.quantities['ymomentum'].vertex_values)
np.save('XstageV.npy', self.domain.quantities['stage'].vertex_values)