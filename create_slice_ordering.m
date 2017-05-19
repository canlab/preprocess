function slice_order_array = create_slice_ordering(slice_order_name, num_slices)

  switch slice_order_name
    case 'interleaved_BU' % "interleaved bottom-up"
      slice_order_array = [1:2:num_slices 2:2:num_slices];
    case 'interleaved_TD' % "interleaved top-down"
      slice_order_array = [num_slices:-2:1, num_slices-1:-2:1];
    case 'ascending'
      slice_order_array = (1:1:num_slices);
    case 'descending'
      slice_order_array = (num_slices:-1:1);
    case 'interleaved_MT'
      for k = 1:num_slices
        slice_order_array(k) = (round((num_slices-k)/2 + (rem((num_slices-k),2) * (num_slices - 1)/2)) + 1);
      end
    otherwise
      error('STOP! Unrecognized image acquisition method');
  end
  
