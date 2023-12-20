

def project = getProject()


path = buildFilePath(PROJECT_BASE_DIR, 'ROIsTif')
mkdirs(path)


for (entry in project.getImageList()) {
    print entry.getImageName()
    
    def baseImageName = entry.getImageName()
    def imageData = entry.readImageData()
    def hierarchy = imageData.getHierarchy()
    
    def server = imageData.getServer()
    println(server.getPath())


    hierarchy.getAnnotationObjects().eachWithIndex{it, x->

      println(" - Working on: "+it)
      def roi = it.getROI()
      def requestROI = RegionRequest.createInstance(server.getPath(), 2.75, roi)
      println(" - checkpoint")
      currentImagePath = buildFilePath('/home/youri/Dropbox/bio-informatica_manuscripts/Article_GLASS-OD/HE-coupes/ROIsTif/', baseImageName + "_"+ x + '.jpg')
      println(" - fn: " + currentImagePath)
      writeImageRegion(server, requestROI, currentImagePath)
        
    }
  
  
}



