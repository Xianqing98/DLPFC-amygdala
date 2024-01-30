function output_clusters = U_clusters_lzr(zs,ps,p_thr)

% p_thr为形成cluster的阈值，可取0.05
% zs为cluster内所有的统计量，ps为相应的p值

output_clusters = [];

cluster_cnt = 0;
p_cnt = 0;
for i = 1:length(zs)
    if(ps(i)<p_thr)
        p_cnt = p_cnt + 1;
    else
        p_cnt = 0;
    end
    if(p_cnt == 2)%create a new cluster
        cluster_cnt = cluster_cnt + 1;
        output_clusters{cluster_cnt}.members = [i-1 i];
        output_clusters{cluster_cnt}.zs = [zs(i-1) zs(i)];
        output_clusters{cluster_cnt}.ps = [ps(i-1) ps(i)];
    elseif(p_cnt > 2)%update the cluster
        output_clusters{cluster_cnt}.members = [output_clusters{cluster_cnt}.members i];
        output_clusters{cluster_cnt}.zs = [output_clusters{cluster_cnt}.zs zs(i)];
        output_clusters{cluster_cnt}.ps = [output_clusters{cluster_cnt}.ps ps(i)];
    end
end

