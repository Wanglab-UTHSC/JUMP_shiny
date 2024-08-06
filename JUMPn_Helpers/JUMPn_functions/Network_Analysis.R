###################################
### Create Node and Edge Frames ###
###################################

Derive_Edge_Weights <- function(Node_Frame, Edge_Frame){
  check_list <- Node_Frame[,1]
  rownames(Node_Frame) <- check_list
  final_edge_frame_2 <- Edge_Frame[((is.element(Edge_Frame[,'from'],check_list)==TRUE) & (is.element(Edge_Frame[,'to'],check_list)==TRUE) ),]
  final_edge_frame_2[,3] <- rep(.1,nrow(final_edge_frame_2))
  final_edge_frame_2[Node_Frame[final_edge_frame_2[,1],3]==Node_Frame[final_edge_frame_2[,2],3],3] <- 2
  final_edge_frame_2 <- cbind(final_edge_frame_2, color='black')
  edge_weights <- unlist(final_edge_frame_2[,3])
  return (list(final_edge_frame_2,edge_weights,Node_Frame, unique(unlist(Node_Frame[,3]))))
}
