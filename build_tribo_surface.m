function z_vals = build_tribo_surface(specify_surface, num_asperities, x_pos, y_pos, zmax, zmin)

    [X,Y] = meshgrid(y_pos, x_pos);
    z_vals = zeros(length(x_pos),length(y_pos));
    rng(0,'twister');
%     rng('shuffle');    
    
    switch specify_surface
        case 1
            % Flat oxide
            z_vals(:,:) = (zmax + zmin)/2.0; %2.5; %0.5*zmin;
        case 2
            % Single node asperities
            z_vals(:,:) = (zmax + zmin)/2.0; %2.5; %0.5*zmin;
            
            num_x_nodes = length(x_pos) - 1;
            num_y_nodes = length(y_pos) - 1;
            
            amount_to_vary = rand((num_x_nodes+1),(num_y_nodes+1)) - 0.5;
            for i = 1:(num_y_nodes+1)
                for j = 1:(num_x_nodes+1)
                    z_vals(j,i) = z_vals(j,i) + amount_to_vary(j,i);
                end
            end
            
        case 3           
            % One paraboloid asperity
            %=====================================================================
            % Define asperity parameters
            %=====================================================================
            z_vals(:,:) = zmin; %2.5; %0.5*zmin;
            
            num_x_nodes = length(x_pos) - 1;
            num_y_nodes = length(y_pos) - 1;
            dx = x_pos(end)/num_x_nodes;
            dy = y_pos(end)/num_y_nodes;
            num_asperities = 1;

            % Define the numnber of nodes in the x-direction that will comprise the
            % asperity
            num_x_nodes_asp = round(0.1 * num_x_nodes); %10;% 7; %5;

            % Define the numnber of nodes in the y-direction that will comprise the
            % asperity    
            num_y_nodes_asp = round(0.1 * num_y_nodes); %10; %7; %5;     


            % This creates a circular asperity inside the bounding box
            asp_radius = 0.9*((num_x_nodes_asp + num_y_nodes_asp)/2.0)*sqrt(dx^2 + dy^2);

            highx = (num_x_nodes+1)-num_x_nodes_asp;
            lowx = num_x_nodes_asp + 1;
            asperity_xs = round(num_x_nodes/2);
            
            highy = (num_y_nodes+1)-num_y_nodes_asp;
            lowy = num_y_nodes_asp + 1;    
            asperity_ys = round(num_y_nodes/2);
            %=====================================================================

            %=====================================================================
            % Create a paraboloid shape that will be used for each asperity
            %=====================================================================
            paraboloid = zeros(2*num_x_nodes_asp + 1, 2*num_y_nodes_asp + 1);
            for j = 1:2*num_x_nodes_asp + 1
                for k = 1:2*num_y_nodes_asp + 1
                    xval = ((j-1)*dx) - (num_x_nodes_asp*dx);
                    yval = ((k-1)*dy) - (num_y_nodes_asp*dy);

                    if sqrt(xval^2 + yval^2) <= asp_radius
                        paraboloid(j,k) = -((xval^2)/(num_x_nodes_asp*dx)^2) - ((yval^2)/(num_y_nodes_asp*dy)^2);
                    end
                end
            end

            hp = max(paraboloid);
            high_p = max(hp);
            lp = min(paraboloid);    
            low_p = min(lp);

            for j = 1:2*num_x_nodes_asp + 1
                for k = 1:2*num_y_nodes_asp + 1
                    xval = ((j-1)*dx) - (num_x_nodes_asp*dx);
                    yval = ((k-1)*dy) - (num_y_nodes_asp*dy);

                    if sqrt(xval^2 + yval^2) <= asp_radius
                        paraboloid(j,k) = paraboloid(j,k) + (-low_p); 
                    end

                end
            end
            %=====================================================================

            %=====================================================================
            % Place each asperity on the oxide surface
            %=====================================================================
            asperity_heights = zmax; %(zmax-zmin).*rand(num_asperities,1) + zmin;

            for i = 1:num_asperities
                idx_x = asperity_xs(i);
                idx_y = asperity_ys(i);

                scale_p = (1.0/(high_p - low_p)) * asperity_heights(i);

                x_asp_low = idx_x - num_x_nodes_asp;
                x_asp_high = idx_x + num_x_nodes_asp;
                xcounter = 1;

                for j = x_asp_low:x_asp_high
                    y_asp_low = idx_y - num_y_nodes_asp;
                    y_asp_high = idx_y + num_y_nodes_asp;
                    ycounter = 1;

                    for k = y_asp_low:y_asp_high
                        z_vals(j,k) = z_vals(j,k) + (scale_p * paraboloid(xcounter,ycounter));
                        ycounter = ycounter + 1;
                    end
                    xcounter = xcounter + 1;
                end

        %         z_vals(idx_x,idx_y) = asperity_heights(i);      
            end
            %=====================================================================            
        case 4
            % Random assortment of paraboloid asperities
            %=====================================================================
            % Define asperity parameters
            %=====================================================================
            num_x_nodes = length(x_pos) - 1;
            num_y_nodes = length(y_pos) - 1;
            dx = x_pos(end)/num_x_nodes;
            dy = y_pos(end)/num_y_nodes;

            % Define the numnber of nodes in the x-direction that will comprise the
            % asperity
            num_x_nodes_asp = round(0.1 * num_x_nodes); %10;% 7; %5;

            % Define the numnber of nodes in the y-direction that will comprise the
            % asperity    
            num_y_nodes_asp = round(0.1 * num_y_nodes); %10; %7; %5;     


            % This creates a circular asperity inside the bounding box
            asp_radius = 0.9*((num_x_nodes_asp + num_y_nodes_asp)/2.0)*sqrt(dx^2 + dy^2);

            highx = (num_x_nodes+1)-num_x_nodes_asp;
            lowx = num_x_nodes_asp + 1;
            asperity_xs = round((highx-lowx).*rand(num_asperities,1) + lowx);
            highy = (num_y_nodes+1)-num_y_nodes_asp;
            lowy = num_y_nodes_asp + 1;    
            asperity_ys = round((highy-lowy).*rand(num_asperities,1) + lowy);
            %=====================================================================

            %=====================================================================
            % Create a paraboloid shape that will be used for each asperity
            %=====================================================================
            paraboloid = zeros(2*num_x_nodes_asp + 1, 2*num_y_nodes_asp + 1);
            for j = 1:2*num_x_nodes_asp + 1
                for k = 1:2*num_y_nodes_asp + 1
                    xval = ((j-1)*dx) - (num_x_nodes_asp*dx);
                    yval = ((k-1)*dy) - (num_y_nodes_asp*dy);

                    if sqrt(xval^2 + yval^2) <= asp_radius
                        paraboloid(j,k) = -((xval^2)/(num_x_nodes_asp*dx)^2) - ((yval^2)/(num_y_nodes_asp*dy)^2);
                    end
                end
            end

            hp = max(paraboloid);
            high_p = max(hp);
            lp = min(paraboloid);    
            low_p = min(lp);

            for j = 1:2*num_x_nodes_asp + 1
                for k = 1:2*num_y_nodes_asp + 1
                    xval = ((j-1)*dx) - (num_x_nodes_asp*dx);
                    yval = ((k-1)*dy) - (num_y_nodes_asp*dy);

                    if sqrt(xval^2 + yval^2) <= asp_radius
                        paraboloid(j,k) = paraboloid(j,k) + (-low_p); 
                    end

                end
            end
            %=====================================================================

            %=====================================================================
            % Place each asperity on the oxide surface
            %=====================================================================
            asperity_heights = (zmax-zmin).*rand(num_asperities,1) + zmin;

            for i = 1:num_asperities
                idx_x = asperity_xs(i);
                idx_y = asperity_ys(i);

                scale_p = (1.0/(high_p - low_p)) * asperity_heights(i);

                x_asp_low = idx_x - num_x_nodes_asp;
                x_asp_high = idx_x + num_x_nodes_asp;
                xcounter = 1;

                for j = x_asp_low:x_asp_high
                    y_asp_low = idx_y - num_y_nodes_asp;
                    y_asp_high = idx_y + num_y_nodes_asp;
                    ycounter = 1;

                    for k = y_asp_low:y_asp_high
                        z_vals(j,k) = z_vals(j,k) + (scale_p * paraboloid(xcounter,ycounter));
                        ycounter = ycounter + 1;
                    end
                    xcounter = xcounter + 1;
                end

        %         z_vals(idx_x,idx_y) = asperity_heights(i);      
            end
            %=====================================================================
            
    end

end