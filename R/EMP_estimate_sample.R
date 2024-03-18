# EMP_estimate_sample <- function(beta_obj,top_number = NULL){
#   pc_data  <- beta_obj$pc_data

#   # 计算近点中心点
#   core_point <- pc_data %>%
#     group_by(Group) %>%
#     summarise(median_PC1_within = median(PC1),median_PC2_within = median(PC2)) %>% dplyr::full_join(.,pc_data,by = 'Group')

#   # 计算远点中心并于近点中心合并
#   pre_data <- pc_data %>%
#     group_by(Group) %>%
#     summarise(median_PC1_without = median(PC1),median_PC2_without = median(PC2)) %>%
#     mutate(Group = recode(Group, "Control" = "Insomnia", "Insomnia" = "Control")) %>% dplyr::full_join(.,core_point,by = 'Group')


#   result <- pre_data %>%
#     mutate(dis_within = c(sqrt( (PC1-median_PC1_within)**2 + (PC2-median_PC2_within)**2 )),
#            dis_without = c(sqrt( (PC1-median_PC1_without)**2 + (PC2-median_PC2_without)**2 )) ) %>%
#     mutate(dis_fold= log(abs(dis_without/dis_within))) %>%
#     mutate(vector1_start = c((PC1 - median_PC1_within)),
#            vector1_end = c((PC2 - median_PC2_within)),
#            vector2_start = c((PC1 - median_PC1_without)),
#            vector2_end = c((PC2 - median_PC2_without))) %>%
#     rowwise() %>%
#     mutate(dis_cosine = cosine(c(vector1_start,vector1_end),c(vector2_start,vector2_end))) %>%
#     mutate(dis_index = dis_fold * dis_cosine) %>%
#     mutate(dis_index = ifelse(dis_fold < 0 | dis_cosine < 0, -abs(dis_index), dis_index)) %>%
#     select(primary,Group,dis_index,dis_fold,dis_cosine,everything())

#   if (is.null(top_number)) {
#     return(result)
#   }else {
#     result %<>% group_by(Group) %>%
#       dplyr::top_n(top_number, dis_index)
#     return(result)
#   }
# }




# 计算向量之间的夹角余弦值
#' Title
#'
#' @param vector1 wait_for_add
#' @param vector2 wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
cosine <- function(vector1,vector2) {
    sum(vector1 * vector2) / (sqrt(sum(vector1^2)) * sqrt(sum(vector2^2)))
}



